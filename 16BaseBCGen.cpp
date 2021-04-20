/*
   Copyright (c) 2021 Ed Harry, Wellcome Sanger Institute, Genome Research Limited

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   */

#define ProgramName "16BaseBCGen"

#include "Common.cpp"

#define ProgramVersion String(PV)

global_function
u32
PackBarCode(u08 *buff)
{
    if (    (buff[1] == '0' && buff[2] == '0') ||
            (buff[4] == '0' && buff[5] == '0') ||
            (buff[7] == '0' && buff[8] == '0') ||
            (buff[10] == '0' && buff[11] == '0')) return(0);

    return( (u32)(((buff[1] - '0') * 10) + (buff[2] - '0')) << 24 |
            (u32)(((buff[4] - '0') * 10) + (buff[5] - '0')) << 16 |
            (u32)(((buff[7] - '0') * 10) + (buff[8] - '0')) << 8  |
            (u32)(((buff[10] - '0') * 10) + (buff[11] - '0')));
}

global_function
void
UnpackBarCode(u32 barcode, u08 *buff)
{
    u08 a = (u08)((barcode >> 24) & ((1 << 8) - 1));
    u08 c = (u08)((barcode >> 16) & ((1 << 8) - 1));
    u08 b = (u08)((barcode >> 8) & ((1 << 8) - 1));
    u08 d = (u08)(barcode & ((1 << 8) - 1));
    stbsp_snprintf((char *)buff, 13, "A%02uC%02uB%02uD%02u", a, c, b, d);
}

global_function
void
Unpack16BaseBarCode(u32 barcode, u08 *buff)
{
    ForLoop(16) *(buff++) = "ATGC"[(barcode >> (2 * (15 - index))) & 3];
}

struct
transfer_buffer_pool
{
    buffer_pool bufferPool;
    wavl_tree *tree;
    memory_arena *arena;
};

global_function
transfer_buffer_pool *
CreateTransferPool(memory_arena *arena, wavl_tree *tree)
{
    transfer_buffer_pool *pool = PushStructP(arena, transfer_buffer_pool);
    pool->bufferPool.pool = ThreadPoolInit(arena, 1);
    pool->tree = tree;
    pool->arena = arena;

    pool->bufferPool.bufferPtr = 0;
    pool->bufferPool.buffers[0] = PushStructP(arena, buffer);
    pool->bufferPool.buffers[0]->buffer = PushArrayP(arena, u08, BufferSize);
    pool->bufferPool.buffers[0]->size = 0;
    pool->bufferPool.buffers[1] = PushStructP(arena, buffer);
    pool->bufferPool.buffers[1]->buffer = PushArrayP(arena, u08, BufferSize);
    pool->bufferPool.buffers[1]->size = 0;

    return(pool);
}

global_function
void
ProcessBuffer(void *in)
{
    transfer_buffer_pool *pool = (transfer_buffer_pool *)in;
    buffer *buffer = pool->bufferPool.buffers[pool->bufferPool.bufferPtr];
    ForLoop64(buffer->size / 4)
    {
        u32 barcode;
        ((u08 *)&barcode)[0] = buffer->buffer[4 * index];
        ((u08 *)&barcode)[1] = buffer->buffer[(4 * index) + 1];
        ((u08 *)&barcode)[2] = buffer->buffer[(4 * index) + 2];
        ((u08 *)&barcode)[3] = buffer->buffer[(4 * index) + 3];

        WavlTreeInsertValue(pool->arena, pool->tree, barcode, 0);
    }
}

global_function
buffer *
GetNextTransferBuffer(transfer_buffer_pool *pool)
{
    FenceIn(ThreadPoolWait(pool->bufferPool.pool));
    buffer *buffer = pool->bufferPool.buffers[pool->bufferPool.bufferPtr];
    pool->bufferPool.bufferPtr = (pool->bufferPool.bufferPtr + 1) & 1;
    ThreadPoolAddTask(pool->bufferPool.pool, ProcessBuffer, pool);
    buffer->size = 0;
    return(buffer);
}

MainArgs
{
    s32 exitCode = EXIT_SUCCESS;
    u08 logError = 0;
    char *logName = (char *)"HaploTag_to_16BaseBCs";

    if (ArgCount > 1 && AreNullTerminatedStringsEqual((u08 *)"--help", (u08 *)ArgBuffer[1])) 
    {
        fprintf(stderr, ProgramName " " ProgramVersion "\nUsage: <fastq format> | " ProgramName " <prefix>? | <fastq format>\n\n");

        fprintf(stderr, "Reads/writes fastq formatted reads from <stdin>/<stdout>.\n");
        fprintf(stderr, "Any read with a BX SAM tag in its comment field will be prepended by 23 bases; a 16-base barcode and 7 joining bases.\n\n");

        fprintf(stderr, "BX tags must be valid haplotag barcodes of the form /^A\\d\\dC\\d\\dB\\d\\dD\\d\\d$/.\n");
        fprintf(stderr, "e.g. '... BX:Z:A01C02B03D04 ...'\n\n");

        fprintf(stderr, "One log file: '%s' will created with an optional '<prefix>_' at the start of the file-name if supplied as an argument.\n", logName);
        fprintf(stderr, "The log file is a map between haplotag and 16-base barcodes.\n");
        fprintf(stderr, "Run 'cut -f 2 HaploTag_to_16BaseBCs | tail -n +2 >16BaseBCs' to extract a list of barcodes suitable for passing as a substitute for a barcode whitelist to other programs.\n\n");

        fprintf(stderr, "Usage example:\n");
        fprintf(stderr, "samtools fastq -@ 16 -nT BX -0 /dev/null -s /dev/null tagged_reads_123.cram | " ProgramName " 123 | bgzip -@ 16 >16BaseBC_reads_123.fq.gz\n");

        goto End;
    }

    char logNameBuffer[256];
    if (ArgCount > 1)
    {
        stbsp_snprintf((char *)logNameBuffer, (s32)sizeof(logNameBuffer), "%s_%s", ArgBuffer[1], logName);
        logName = (char *)logNameBuffer;
    }

    s32 log;
    if ((log = open((const char *)logName, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR)) > 0)
    {
        memory_arena workingSet;
        CreateMemoryArena(workingSet, MegaByte(512));
        
        wavl_tree *barcodeTree = InitialiseWavlTree(&workingSet);
        transfer_buffer_pool *transferBufferPool = CreateTransferPool(&workingSet, barcodeTree);
        
        buffer_pool *readPool = CreatePool(&workingSet);
#ifdef DEBUG
        readPool->handle = open("test_in", O_RDONLY);
#else
        readPool->handle = STDIN_FILENO;
#endif     
        buffer_pool *writePool = CreatePool(&workingSet);
        writePool->handle = STDOUT_FILENO;

        char printNBuffers[2][32] = {{0}};
        u08 printNBufferPtr = 0;

        enum modes {null, readTag1, readTag2, readTag3, readTag4, readTag5, readingData, gotBC, write1, doneWrite1, skip1, write2, skip2, skip3, done};
        modes mode = null;
        u08 bcBuffer[13];
        u08 bcPtr = 0;
        u32 barcode = 0;

        buffer *writeBuffer = GetNextBuffer_Write(writePool);
        buffer *readBuffer = GetNextBuffer_Read(readPool);
        buffer *transferBuffer = GetNextTransferBuffer(transferBufferPool);
        u64 totalReads = 0;
        u64 bcAdded = 0;
        do
        {
            readBuffer = GetNextBuffer_Read(readPool);

            for (   u64 bufferIndex = 0;
                    bufferIndex < readBuffer->size;
                    ++bufferIndex )
            {
                if (Global_Write_Error)
                {
                    PrintError("Error writing");
                    exitCode = EXIT_FAILURE;
                    goto End;
                }

                u08 character = readBuffer->buffer[bufferIndex];

                if (mode == null && character < 33) mode = readTag1;
                else if (mode == readTag1) mode = character == 'B' ? readTag2 : null;
                else if (mode == readTag2) mode = character == 'X' ? readTag3 : null;
                else if (mode == readTag3) mode = character == ':' ? readTag4 : null;
                else if (mode == readTag4) mode = character == 'Z' ? readTag5 : null;
                else if (mode == readTag5) 
                {
                    mode = character == ':' ? readingData : null;
                    bcPtr = 0;
                }
                else if (mode == readingData)
                {
                    if (character < 33) mode = bcPtr == (sizeof(bcBuffer) - 1) ? gotBC : null;
                    else
                    {
                        bcBuffer[bcPtr++] = character;
                        if (bcPtr == sizeof(bcBuffer)) mode = null;
                    }
                }

                if (character == '\n')
                {
                    if (mode == done) 
                    {
                        mode = null;
                        ++totalReads;
                    }
                    else if (mode == skip3) mode = done;
                    else if (mode == skip2) mode = skip3; 
                    else if (mode == skip1) mode = write2;
                    else if (mode == doneWrite1) mode = skip1;
                    else if (mode == gotBC) mode = write1;
                    else mode = skip2; 
                }
                else if (mode == write1)
                {
                    if ((BufferSize - writeBuffer->size - 1) < 23) writeBuffer = GetNextBuffer_Write(writePool);
                    
                    barcode = PackBarCode(bcBuffer);
                    u08 tenx[16];

                    if (barcode)
                    {
                        transferBuffer->buffer[transferBuffer->size++] = ((u08 *)&barcode)[0];
                        transferBuffer->buffer[transferBuffer->size++] = ((u08 *)&barcode)[1];
                        transferBuffer->buffer[transferBuffer->size++] = ((u08 *)&barcode)[2];
                        transferBuffer->buffer[transferBuffer->size++] = ((u08 *)&barcode)[3];
                        if (transferBuffer->size == BufferSize) transferBuffer = GetNextTransferBuffer(transferBufferPool);
                        
                        Unpack16BaseBarCode(barcode, tenx);
                    }

                    ForLoop(16) writeBuffer->buffer[writeBuffer->size++] = barcode ? tenx[index] : 'N';
                    ForLoop(7) writeBuffer->buffer[writeBuffer->size++] = barcode ? 'A' : 'N';
                    mode = doneWrite1;
                }
                else if (mode == write2)
                {
                    if ((BufferSize - writeBuffer->size - 1) < 23) writeBuffer = GetNextBuffer_Write(writePool);
                    
                    ForLoop(23) writeBuffer->buffer[writeBuffer->size++] = barcode ? 'J' : '#';
                    ++bcAdded;
                    
                    mode = done;
                }

                writeBuffer->buffer[writeBuffer->size++] = character;
                if (BufferSize == writeBuffer->size) writeBuffer = GetNextBuffer_Write(writePool);

#define Log2_Print_Interval 14
                if ((totalReads && !(totalReads & ((1 << Log2_Print_Interval) - 1))) || (bcAdded && !(bcAdded & ((1 << Log2_Print_Interval) - 1))))
                {
                    u08 currPtr = printNBufferPtr;
                    u08 otherPtr = (currPtr + 1) & 1;
                    stbsp_snprintf(printNBuffers[currPtr], sizeof(printNBuffers[currPtr]), "%$" PRIu64 " / %$" PRIu64, totalReads, bcAdded);

                    if (strcmp(printNBuffers[currPtr], printNBuffers[otherPtr]))
                    {
                        PrintStatus("%s reads processed / barcodes added", printNBuffers[currPtr]);
                    }

                    printNBufferPtr = otherPtr;
                }
            }
        } while (readBuffer->size);

        GetNextBuffer_Write(writePool);
        if (Global_Write_Error)
        {
            PrintError("Error writing");
            exitCode = EXIT_FAILURE;
        }      

        {
            GetNextTransferBuffer(transferBufferPool);
            GetNextTransferBuffer(transferBufferPool);
            WavlTreeFreeze_LowToHigh(barcodeTree);

            char *header = (char *)"HaploTag\t16 Base BC\n";
            if (WriteToLogFile(log, header, strlen(header)))
            {
                logError = 1;
                goto End;
            }

            TraverseLinkedList(WavlTreeGetBottom(barcodeTree)->next, wavl_node)
            {
                u08 lineBuffer[30];
                UnpackBarCode(node->value, lineBuffer);
                lineBuffer[12] = '\t';
                Unpack16BaseBarCode(node->value, lineBuffer + 13);
                lineBuffer[29] = '\n';

                if (WriteToLogFile(log, (char *)lineBuffer, 30))
                {
                    logError = 1;
                    goto End;
                }
            }  
        }

        GetNextBuffer_Write(writePool);
        if (Global_Write_Error)
        {
            PrintError("Error writing");
            exitCode = EXIT_FAILURE;
        }
    }
    else
    {
        PrintError("Error opening log file");
        exitCode = EXIT_FAILURE;
    }

End:
    if (logError)
    {
        PrintError("Error writing log file");
        exitCode = EXIT_FAILURE;
    }

    return(exitCode);
}
