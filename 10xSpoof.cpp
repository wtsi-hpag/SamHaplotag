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

#define ProgramName "10xSpoof"

#include "Common.cpp"

#define ProgramVersion String(PV)

#include "10x.cpp"

struct
barcode
{
    u32 barcode;
    u32 index;
};

struct
barcode_hash_table_node
{
    barcode *barcode;
    barcode_hash_table_node *next;
};

struct
barcode_hash_table
{
    barcode_hash_table_node **table;
    u32 size;
    u32 pad;
};

global_function
barcode_hash_table *
CreateBarCodeHashTable(memory_arena *arena, u32 nEntries)
{
    u32 size = NextPrime((u32)((f32)nEntries * 1.3f));

    barcode_hash_table *table = PushStructP(arena, barcode_hash_table);
    table->size = size;
    table->table = PushArrayP(arena, barcode_hash_table_node *, size);
    memset(table->table, 0, size * sizeof(barcode_hash_table_node *));

    return(table);
}

global_function
void
AddBarCodeToHashTable(barcode_hash_table *table, memory_arena *arena, u32 code, u32 index)
{
#define BarCodeHashTableSeed 0xf40b503a1896576e
    u32 hash = FastHash32(&code, sizeof(code), BarCodeHashTableSeed) % table->size;
    barcode_hash_table_node *node = table->table[hash];
    barcode_hash_table_node *prevNode = 0;

    while (node)
    {
        prevNode = node;
        node = node->next;
    }

    barcode_hash_table_node *newNode = PushStructP(arena, barcode_hash_table_node);
    newNode->barcode = PushStructP(arena, barcode);
    newNode->barcode->barcode = code;
    newNode->barcode->index = index;
    newNode->next = 0;

    (prevNode ? prevNode->next : table->table[hash]) = newNode;
}

global_function
u32
GetBarCodeIndexFromHashTable(barcode_hash_table *table, u32 barcode)
{
    u32 hash = FastHash32(&barcode, sizeof(barcode), BarCodeHashTableSeed) % table->size;
    barcode_hash_table_node *node = table->table[hash];
    u32 result = 0;

    while (node)
    {
        if (node->barcode->barcode == barcode)
        {
            result = node->barcode->index;
            break;
        }
        node = node->next;
    }

    return(result);
}

global_function
u32
PackBarCode(u08 *buff)
{
    return( (u32)(((buff[1] - '0') * 10) + (buff[2] - '0')) << 24 |
            (u32)(((buff[4] - '0') * 10) + (buff[5] - '0')) << 16 |
            (u32)(((buff[7] - '0') * 10) + (buff[8] - '0')) << 8  |
            (u32)(((buff[10] - '0') * 10) + (buff[11] - '0')));
}

MainArgs
{
    s32 exitCode = EXIT_SUCCESS;
    u08 logError = 0;
    char *logName = (char *)"10xSpoof_HaploTag_to_10x";
    
    if (ArgCount > 1 && AreNullTerminatedStringsEqual((u08 *)"--help", (u08 *)ArgBuffer[1])) 
    {
        fprintf(stderr, ProgramName " " ProgramVersion "\nUsage: <fastq format> | " ProgramName " <clear barcode log> <prefix>? | <fastq format>\n\n");
        
        fprintf(stderr, "Reads/writes fastq formatted reads from <stdin>/<stdout>.\n");
        fprintf(stderr, "Any read with a BX SAM tag in its comment field will be prepended by 23 bases; a 16-base valid 10x barcode and 7 joining bases.\n\n");
        
        fprintf(stderr, "BX tags must be valid haplotag barcodes of the form /^A\\d\\dC\\d\\dB\\d\\dD\\d\\d$/.\n");
        fprintf(stderr, "e.g. '... BX:Z:A01C02B03D04 ...'\n\n");
       
        fprintf(stderr, "The one required argument <clear barcode log> is a 3-column, tab-delimited text file with one header line; with the columns being: haplotag barcode, clear-count and correct-count.\n");
        fprintf(stderr, "Such a log file will be created by running 'SamHaplotag'.\n\n");

        fprintf(stderr, "One log file: '%s' will created with an optional '<prefix>_' at the start of the file-name if supplied as a second argument.\n", logName);
        fprintf(stderr, "The log file is a map between haplotag and 10x barcodes.\n\n");

        fprintf(stderr, "Usage example:\n");
        fprintf(stderr, "samtools fastq -@ 16 -nT BX -0 /dev/null -s /dev/null tagged_reads_123.cram | " ProgramName " 123_SamHaplotag_Clear_BC 123 | bgzip -@ 16 >10x_spoofed_reads_123.fq.gz\n");
        
        goto End;
    }

    if (ArgCount < 2)
    {
        PrintError("Clear Barcode log required");
        exitCode = EXIT_FAILURE;
    }
    else
    {
        char logNameBuffer[256];
        if (ArgCount > 2)
        {
            stbsp_snprintf((char *)logNameBuffer, (s32)sizeof(logNameBuffer), "%s_%s", ArgBuffer[2], logName);
            logName = (char *)logNameBuffer;
        }

        s32 log;
        if ((log = open((const char *)logName, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR)) > 0)
        {
            PrintStatus("Clear Barcode log: %s", ArgBuffer[1]);

            memory_arena workingSet;
            CreateMemoryArena(workingSet, MegaByte(512));

            buffer_pool *readPool = CreatePool(&workingSet);
            readPool->handle = open(ArgBuffer[1], O_RDONLY);

            if (readPool->handle > 0)
            {
                barcode_hash_table *barcodeHashTable;
                wavl_tree *barcodeTree = InitialiseWavlTree(&workingSet);
                {
                    u32 barcode = 0;
                    u32 count = 0;

                    enum logState {head, bc, n1, n2};
                    logState state = head;

                    u08 buff[16];
                    u08 buffPtr = 0;

                    buffer *readBuffer = GetNextBuffer_Read(readPool);
                    u32 nBC = 0;
                    do
                    {
                        readBuffer = GetNextBuffer_Read(readPool);

                        for (   u64 bufferIndex = 0;
                                bufferIndex < readBuffer->size;
                                ++bufferIndex )
                        {
                            u08 character = readBuffer->buffer[bufferIndex];

                            switch (state)
                            {
                                case head:
                                    if (character == '\n') state = bc;
                                    break;

                                case bc:
                                    if (character == '\t')
                                    {
                                        barcode = PackBarCode(buff);
                                        state = n1;
                                        buffPtr = 0;
                                    }
                                    else buff[buffPtr++] = character;
                                    break;

                                case n1:
                                    if (character == '\t')
                                    {
                                        count = StringToInt(buff + buffPtr, buffPtr);
                                        state = n2;
                                        buffPtr = 0;
                                    }
                                    else buff[buffPtr++] = character;
                                    break;

                                case n2:
                                    if (character == '\n')
                                    {
                                        count += StringToInt(buff + buffPtr, buffPtr);
                                        state = bc;
                                        buffPtr = 0;

                                        WavlTreeInsertValue(&workingSet, barcodeTree, count, barcode);
                                        ++nBC;
                                    }
                                    else buff[buffPtr++] = character;
                            }
                        }
                    } while (readBuffer->size);

                    WavlTreeFreeze_HighToLow(barcodeTree);

                    PrintStatus("Barcode count: %u", nBC);
                    barcodeHashTable = CreateBarCodeHashTable(&workingSet, nBC);

                    if (nBC > ArrayCount(TenX_BarCodes)) PrintWarning("Barcode count > 10x count, %u barcodes will be discarded!", nBC - ArrayCount(TenX_BarCodes));

                    char *header = (char *)"HaploTag\t10x\n";
                    if (WriteToLogFile(log, header, strlen(header)))
                    {
                        logError = 1;
                        goto End;
                    }

                    u32 index = 1;
                    TraverseLinkedList(WavlTreeGetTop(barcodeTree), wavl_node)
                    {
                        if (!node->value) break;
                        TraverseLinkedList2(node->head, barcode_ll_node)
                        {
                            u08 lineBuffer[30];
                            u08 a = (u08)((node2->barcode >> 24) & ((1 << 8) - 1));
                            u08 c = (u08)((node2->barcode >> 16) & ((1 << 8) - 1));
                            u08 b = (u08)((node2->barcode >> 8) & ((1 << 8) - 1));
                            u08 d = (u08)(node2->barcode & ((1 << 8) - 1));
                            stbsp_snprintf((char *)lineBuffer, sizeof(lineBuffer), "A%02uC%02uB%02uD%02u\t", a, c, b, d);
                            UnPackTenX(index - 1, lineBuffer + 13);
                            lineBuffer[29] = '\n';

                            if (WriteToLogFile(log, (char *)lineBuffer, 30))
                            {
                                logError = 1;
                                goto End;
                            }

                            AddBarCodeToHashTable(barcodeHashTable, &workingSet, node2->barcode, index++);

                            if (index > ArrayCount(TenX_BarCodes)) break;
                        }
                        if (index > ArrayCount(TenX_BarCodes)) break;
                    }
                }
                {
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
                    u32 bcIndex = 0;

                    buffer *writeBuffer = GetNextBuffer_Write(writePool);
                    buffer *readBuffer = GetNextBuffer_Read(readPool);
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
                                bcIndex = GetBarCodeIndexFromHashTable(barcodeHashTable, PackBarCode(bcBuffer));
                                u08 tenx[16];
                                if (bcIndex) UnPackTenX(bcIndex - 1, tenx);
                                ForLoop(16) writeBuffer->buffer[writeBuffer->size++] = bcIndex ? tenx[index] : 'N';
                                ForLoop(7) writeBuffer->buffer[writeBuffer->size++] = bcIndex ? 'A' : 'N';
                                mode = doneWrite1;
                            }
                            else if (mode == write2)
                            {
                                if ((BufferSize - writeBuffer->size - 1) < 23) writeBuffer = GetNextBuffer_Write(writePool);
                                ForLoop(23) writeBuffer->buffer[writeBuffer->size++] = bcIndex ? 'J' : '#';
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
                    GetNextBuffer_Write(writePool);
                    if (Global_Write_Error)
                    {
                        PrintError("Error writing");
                        exitCode = EXIT_FAILURE;
                    }
                } 
            }
            else
            {
                PrintError("Error opening log file '%s'", ArgBuffer[1]);
                exitCode = EXIT_FAILURE;
            }
        }
        else
        {
            PrintError("Error opening log file");
            exitCode = EXIT_FAILURE;
        }
    }

End:
    if (logError)
    {
        PrintError("Error writing log file");
        exitCode = EXIT_FAILURE;
    }
    
    return(exitCode);
}
