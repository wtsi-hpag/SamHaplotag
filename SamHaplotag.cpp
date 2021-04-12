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

#define ProgramName "SamHaplotag"
#define ProgramVersion "0.0.2"

#include "Common.cpp"
#include "BC.cpp"

global_function
u08
Comp(u08 base)
{
    u08 table[] = {'T', 0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0, 0, 0, 0, 0, 'A'};
    return(table[base-65]);
}

struct
string_hash_table_node
{
    u08 *string;
    string_hash_table_node *next;
};

struct
string_hash_table
{
    string_hash_table_node **table;
    u32 size;
    u32 pad;
};

global_function
string_hash_table *
CreateStringHashTable(memory_arena *arena)
{
    u32 size = 131;

    string_hash_table *table = PushStructP(arena, string_hash_table);
    table->size = size;
    table->table = PushArrayP(arena, string_hash_table_node *, size);
    memset(table->table, 0, size * sizeof(string_hash_table_node *));

    return(table);
}

global_function
void
AddStringToHashTable(string_hash_table *table, memory_arena *arena, u08 *string)
{
#define StringHashTableSeed 0xf30b603a1896576e
    u32 hash = FastHash32(string, StringLength(string), StringHashTableSeed) % table->size;
    string_hash_table_node *node = table->table[hash];
    string_hash_table_node *prevNode = 0;

    while (node)
    {
        prevNode = node;
        node = node->next;
    }

    string_hash_table_node *newNode = PushStructP(arena, string_hash_table_node);
    newNode->string = string;
    newNode->next = 0;
    (prevNode ? prevNode->next : table->table[hash]) = newNode;
}

global_function
u08
IsStringInHashTable(string_hash_table *table, memory_arena *arena, u08 *string)
{
    u32 hash = FastHash32(string, StringLength(string), StringHashTableSeed) % table->size;
    string_hash_table_node *node = table->table[hash];
    u08 result = 0;

    while (node)
    {
        if (AreNullTerminatedStringsEqual(string, node->string))
        {
            result = 1;
            break;
        }
        node = node->next;
    }

    return(result);
}

struct
barcode
{
    u32 barcode;
    u32 correct;
    u32 corrected;
    u32 unclear;
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
CreateBarCodeHashTable(memory_arena *arena)
{
    u32 size = 0x14CCCF;

    barcode_hash_table *table = PushStructP(arena, barcode_hash_table);
    table->size = size;
    table->table = PushArrayP(arena, barcode_hash_table_node *, size);
    memset(table->table, 0, size * sizeof(barcode_hash_table_node *));

    return(table);
}

global_function
barcode *
GetBarCodeFromHashTable(barcode_hash_table *table, memory_arena *arena, u32 code)
{
#define BarCodeHashTableSeed 0xf30b503a1896576e
    u32 hash = FastHash32(&code, sizeof(code), BarCodeHashTableSeed) % table->size;
    barcode_hash_table_node *node = table->table[hash];
    barcode_hash_table_node *prevNode = 0;
    barcode *result = 0;

    while (node)
    {
        prevNode = node;
        if (code == node->barcode->barcode)
        {
            result = node->barcode;
            break;
        }
        node = node->next;
    }

    if (!result)
    {
        barcode_hash_table_node *newNode = PushStructP(arena, barcode_hash_table_node);
        newNode->barcode = PushStructP(arena, barcode);
        newNode->next = 0;
        newNode->barcode->barcode = code;
        newNode->barcode->correct = 0;
        newNode->barcode->corrected = 0;
        newNode->barcode->unclear = 0;
        (prevNode ? prevNode->next : table->table[hash]) = newNode;

        result = newNode->barcode;
    }

    return(result);
}

struct
transfer_buffer_pool
{
    buffer_pool bufferPool;
    barcode_hash_table *table;
    wavl_tree *tree;
    memory_arena *arena;
};

global_function
transfer_buffer_pool *
CreateTransferPool(memory_arena *arena, barcode_hash_table *table, wavl_tree *tree)
{
    transfer_buffer_pool *pool = PushStructP(arena, transfer_buffer_pool);
    pool->bufferPool.pool = ThreadPoolInit(arena, 1);
    pool->table = table;
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
        u08 a = buffer->buffer[4 * index];
        u08 b = buffer->buffer[(4 * index) + 1];
        u08 c = buffer->buffer[(4 * index) + 2];
        u08 d = buffer->buffer[(4 * index) + 3];

        u32 barcodePack = (((u32)(a & 127)) << 24) | (((u32)(c & 127)) << 16) | (((u32)(b & 127)) << 8) | ((u32)(d & 127));
        barcode *barcode = GetBarCodeFromHashTable(pool->table, pool->arena, barcodePack);
        WavlTreeInsertValue(pool->arena, pool->tree, barcodePack + 1, 0);

        ++((a && b && c && d) ? (((a | b | c | d) & 128) ? barcode->corrected : barcode->correct) : barcode->unclear);
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

    char *missingTagsLogName = (char *)"SamHaplotag_Missing_BC_QT_tags";
    char *clearBCLogName = (char *)"SamHaplotag_Clear_BC";
    char *unclearBCLogName = (char *)"SamHaplotag_UnClear_BC";
    
    if (ArgCount > 1 && AreNullTerminatedStringsEqual((u08 *)"--help", (u08 *)ArgBuffer[1])) 
    {
        fprintf(stderr, ProgramName " " ProgramVersion "\nUsage: <sam format> | " ProgramName " <prefix>? | <sam format>\n\n");
        
        fprintf(stderr, "Reads/writes SAM formatted reads from <stdin>/<stdout>.\n");
        fprintf(stderr, "Any reads flagged as <read1> with both BC and QT tags will have additional haplotag RX, QX and BX tags added.\n\n");
        
        fprintf(stderr, "BC tags must be of the form /^[ATGCN]{13}\\-[ATGCN]{13}$/ and QT tags of the form /^[!-~]{13}\\w[!-~]{13}$/.\n");
        fprintf(stderr, "e.g. '... BC:Z:NGGTACATGAGAC-NTATCGGCCTTCA\tQT:Z:!FFFFFFFFFFFF !,,,F,FFF:F:F ...'\n\n");
       
        fprintf(stderr, "Three log files: '%s', '%s' and '%s' are created with an optional '<prefix>_' at the start of each file-name if supplied as an argument.\n\n", clearBCLogName, unclearBCLogName, missingTagsLogName);

        fprintf(stderr, "Usage example:\n");
        fprintf(stderr, "samtools view -h@ 16 -F 0xF00 reads_123.cram | " ProgramName " 123 | samtools view -@ 16 -o tagged_reads_123.cram\n");
        
        goto End;
    }
    
    char logNameBuffer[256];
    if (ArgCount > 1)
    {
        s32 ptr1 = 1 + stbsp_snprintf((char *)logNameBuffer, (s32)sizeof(logNameBuffer), "%s_%s", ArgBuffer[1], missingTagsLogName);
        missingTagsLogName = (char *)logNameBuffer;

        s32 ptr2 = 1 + stbsp_snprintf((char *)logNameBuffer + ptr1, (s32)sizeof(logNameBuffer) - ptr1, "%s_%s", ArgBuffer[1], clearBCLogName);
        clearBCLogName = (char *)logNameBuffer + ptr1;

        stbsp_snprintf((char *)logNameBuffer + ptr1 + ptr2, (s32)sizeof(logNameBuffer) - ptr1 - ptr2, "%s_%s", ArgBuffer[1], unclearBCLogName);
        unclearBCLogName = (char *)logNameBuffer + ptr1 + ptr2;
    }

    PrintStatus("Starting...");
    
    s32 missingTagsLog, clearBCLog, unclearBCLog;
    if (    (missingTagsLog = open((const char *)missingTagsLogName, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR)) > 0 &&
            (clearBCLog = open((const char *)clearBCLogName, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR)) > 0 &&
            (unclearBCLog = open((const char *)unclearBCLogName, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR)) > 0
       )
    {
        char *header = (char *)"Read\n";
        if (WriteToLogFile(missingTagsLog, header, strlen(header)))
        {
            logError = 1;
            goto End;
        }

        memory_arena workingSet;
        CreateMemoryArena(workingSet, MegaByte(512));

        barcode_hash_table *barcodeHashTable = CreateBarCodeHashTable(&workingSet);
        wavl_tree *barcodeTree = InitialiseWavlTree(&workingSet);

        transfer_buffer_pool *transferBufferPool = CreateTransferPool(&workingSet, barcodeHashTable, barcodeTree);

        buffer_pool *readPool = CreatePool(&workingSet);
#ifdef DEBUG
        readPool->handle = open("test_in", O_RDONLY);
#else
        readPool->handle = STDIN_FILENO;
#endif        
        buffer_pool *writePool = CreatePool(&workingSet);
        writePool->handle = STDOUT_FILENO;

        char printNBuffers[2][16] = {{0}};
        u08 printNBufferPtr = 0;

        u08 headerMode = 1;
        u08 atEnd = 1;
        u64 total = 0;

        u08 nameBuffer[64];
        u08 namePtr = 0;

        u32 flags = 0;

        u08 BCBuffer[27];
        u08 QTBuffer[27];
        u08 flagBuffer[5];
        u08 tagPtr = 0;
        enum tagStat {null, readTag1, readTag2, readTag3, readTag4, readTag5, readingData, done};
        tagStat BC = null;
        tagStat QT = null;
        tagStat FL = null;

        u08 IDLine[64];
        tagStat PG = null;
        tagStat ID = null;
        string_hash_table *ids = CreateStringHashTable(&workingSet);
        u08 *lastID = 0;
        //@PG     ID:samtools.21  PN:samtools     PP:samtools.20  VN:1.12 CL:samtools view -H PycharmProjects/HiLine/TestData/valid.cram

        buffer *readBuffer = GetNextBuffer_Read(readPool);
        buffer *writeBuffer = GetNextBuffer_Write(writePool);
        buffer *transferBuffer = GetNextTransferBuffer(transferBufferPool);
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

                if (headerMode && atEnd) 
                {
                    headerMode = character == '@';
                    if (!headerMode)
                    {
                        tagPtr = 0;
                        
                        u08 idBuff[64];
                        stbsp_snprintf((char *)idBuff, sizeof(idBuff), "%s", ProgramName);
                        u32 c = 0;
                        while (IsStringInHashTable(ids, &workingSet, idBuff)) stbsp_snprintf((char *)idBuff, sizeof(idBuff), "%s.%u", ProgramName, ++c);
                        
                        u08 pgLine[512];
                        u32 n = lastID ?    (u32)stbsp_snprintf((char *)pgLine, sizeof(pgLine), "@PG\tID:%s\tPN:%s\tPP:%s\tVN:%s\tCL:", idBuff, ProgramName, (char *)lastID, ProgramVersion) :
                                            (u32)stbsp_snprintf((char *)pgLine, sizeof(pgLine), "@PG\tID:%s\tPN:%s\tVN:%s\tCL:", idBuff, ProgramName, ProgramVersion);
                        
                        ForLoop(ArgCount) n += (u32)stbsp_snprintf((char *)pgLine + n, sizeof(pgLine) - n, "%s ", ArgBuffer[index]);
                        pgLine[n - 1] = '\n';

                        if ((BufferSize - writeBuffer->size - 1) < n) writeBuffer = GetNextBuffer_Write(writePool);
                        ForLoop(n) writeBuffer->buffer[writeBuffer->size++] = pgLine[index];
                    }
                }
                atEnd = character == '\n';

                if (headerMode)
                {
                    if (character == '@') PG = readTag1;
                    else if (PG == readTag1) PG = character == 'P' ? readTag2 : null;
                    else if (PG == readTag2) PG = character == 'G' ? readTag3 : null;
                    else if (PG == readTag3) PG = character == '\t' ? readTag4 : null;

                    if (PG == readTag4 && character == '\t')
                    {
                        ID = readTag1;
                        PG = null;
                    }
                    else if (ID == readTag1) ID = character == 'I' ? readTag2 : null;
                    else if (ID == readTag2) ID = character == 'D' ? readTag3 : null;
                    else if (ID == readTag3) 
                    {
                        ID = character == ':' ? readingData : null;
                        tagPtr = 0;
                    }
                    else if (ID == readingData)
                    {
                        IDLine[tagPtr++] = character;
                        if (character < 33)
                        {
                            ID = done;
                            IDLine[tagPtr-1] = 0;
                            u08 *str = PushArray(workingSet, u08, tagPtr);
                            ForLoop((u32)tagPtr) str[index] = IDLine[index];
                            AddStringToHashTable(ids, &workingSet, str);
                            lastID = str;
                        }
                    }
                }
                else
                {
                    if (namePtr < (sizeof(nameBuffer) - 1)) 
                    {
                        if (character == '\t')
                        {
                            nameBuffer[namePtr] = 0;
                            namePtr = sizeof(nameBuffer);
                        }
                        else
                        {
                            nameBuffer[namePtr++] = character;
                            if (namePtr == (sizeof(nameBuffer) - 1)) nameBuffer[namePtr] = 0;
                        }
                    }

                    if (!atEnd)
                    {
                        if (FL == null && character == '\t') FL = readingData;
                        else if (FL == readingData)
                        {
                            if (character != '\t')
                            {
                                flagBuffer[tagPtr++] = character;
                            }
                            else
                            {
                                flags = StringToInt(flagBuffer + tagPtr, (u32)tagPtr);

                                tagPtr = 0;
                                FL = done;
                            }
                        }

                        if (BC == null && character == '\t') BC = readTag1;
                        else if (BC == readTag1) BC = character == 'B' ? readTag2 : null;
                        else if (BC == readTag2) BC = character == 'C' ? readTag3 : null;
                        else if (BC == readTag3) BC = character == ':' ? readTag4 : null;
                        else if (BC == readTag4) BC = character == 'Z' ? readTag5 : null;
                        else if (BC == readTag5) BC = character == ':' ? readingData : null;
                        else if (BC == readingData)
                        {
                            BCBuffer[tagPtr++] = character;
                            if (tagPtr == sizeof(BCBuffer))
                            {
                                tagPtr = 0;
                                BC = done;
                            }
                        }

                        if (QT == null && character == '\t') QT = readTag1;
                        else if (QT == readTag1) QT = character == 'Q' ? readTag2 : null;
                        else if (QT == readTag2) QT = character == 'T' ? readTag3 : null;
                        else if (QT == readTag3) QT = character == ':' ? readTag4 : null;
                        else if (QT == readTag4) QT = character == 'Z' ? readTag5 : null;
                        else if (QT == readTag5) QT = character == ':' ? readingData : null;
                        else if (QT == readingData)
                        {
                            QTBuffer[tagPtr++] = character;
                            if (tagPtr == sizeof(QTBuffer))
                            {
                                tagPtr = 0;
                                QT = done;
                            }
                        }
                    }
                    else
                    {
                        if (flags & 64)
                        {
                            if (BC == done && QT == done)
                            {
                                u32 totalNewSpace = (2 * (6 + 27)) + 6 + 12;
                                if ((BufferSize - writeBuffer->size - 1) < totalNewSpace) writeBuffer = GetNextBuffer_Write(writePool);

                                u08 revComBuffer[13];
                                ForLoop(13) revComBuffer[index] = Comp(BCBuffer[sizeof(BCBuffer) - index - 1]);

                                // RX
                                writeBuffer->buffer[writeBuffer->size++] = '\t';
                                writeBuffer->buffer[writeBuffer->size++] = 'R';
                                writeBuffer->buffer[writeBuffer->size++] = 'X';
                                writeBuffer->buffer[writeBuffer->size++] = ':';
                                writeBuffer->buffer[writeBuffer->size++] = 'Z';
                                writeBuffer->buffer[writeBuffer->size++] = ':';
                                ForLoop(13) writeBuffer->buffer[writeBuffer->size++] = BCBuffer[index];
                                writeBuffer->buffer[writeBuffer->size++] = '+';
                                ForLoop(13) writeBuffer->buffer[writeBuffer->size++] = revComBuffer[index];

                                // QX
                                writeBuffer->buffer[writeBuffer->size++] = '\t';
                                writeBuffer->buffer[writeBuffer->size++] = 'Q';
                                writeBuffer->buffer[writeBuffer->size++] = 'X';
                                writeBuffer->buffer[writeBuffer->size++] = ':';
                                writeBuffer->buffer[writeBuffer->size++] = 'Z';
                                writeBuffer->buffer[writeBuffer->size++] = ':';
                                ForLoop(13) writeBuffer->buffer[writeBuffer->size++] = QTBuffer[index];
                                writeBuffer->buffer[writeBuffer->size++] = '+';
                                ForLoop(13) writeBuffer->buffer[writeBuffer->size++] = QTBuffer[sizeof(QTBuffer) - index - 1];

                                // BX
                                u08 BXBuffer[13];

                                u08 a = GetBC_A(BCBuffer + 7);
                                u08 b = GetBC_B(revComBuffer + 7);
                                u08 c = GetBC_C(BCBuffer);
                                u08 d = GetBC_D(revComBuffer);

                                transferBuffer->buffer[transferBuffer->size++] = a;
                                transferBuffer->buffer[transferBuffer->size++] = b;
                                transferBuffer->buffer[transferBuffer->size++] = c;
                                transferBuffer->buffer[transferBuffer->size++] = d;
                                if (transferBuffer->size == BufferSize) transferBuffer = GetNextTransferBuffer(transferBufferPool);

                                stbsp_snprintf((char *)BXBuffer, 13, "A%02uC%02uB%02uD%02u", a&127, c&127, b&127, d&127);
                                writeBuffer->buffer[writeBuffer->size++] = '\t';
                                writeBuffer->buffer[writeBuffer->size++] = 'B';
                                writeBuffer->buffer[writeBuffer->size++] = 'X';
                                writeBuffer->buffer[writeBuffer->size++] = ':';
                                writeBuffer->buffer[writeBuffer->size++] = 'Z';
                                writeBuffer->buffer[writeBuffer->size++] = ':';
                                ForLoop(12) writeBuffer->buffer[writeBuffer->size++] = BXBuffer[index];
                            }
                            else
                            {
                                PrintWarning("Read %s has no %s tag%s", nameBuffer, (BC != done && QT != done) ? "BC/QT" : (BC != done ? "BC" : "QT"), (BC != done && QT != done) ? "s" : "");

                                if (WriteToLogFile(missingTagsLog, nameBuffer, strlen((char *)nameBuffer)) || WriteToLogFile(missingTagsLog, (void *)"\n", 1))
                                {
                                    logError = 1;
                                    goto End;
                                }
                            }
                        }
                        
                        FL = BC = QT = null;
                        tagPtr = namePtr = 0;
                    }
                }

                writeBuffer->buffer[writeBuffer->size++] = character;
                if (BufferSize == writeBuffer->size) writeBuffer = GetNextBuffer_Write(writePool);

                if (!headerMode && atEnd)
                {
#define Log2_Print_Interval 14
                    if (!(++total & ((1 << Log2_Print_Interval) - 1)))
                    {
                        u08 currPtr = printNBufferPtr;
                        u08 otherPtr = (currPtr + 1) & 1;
                        stbsp_snprintf(printNBuffers[currPtr], sizeof(printNBuffers[currPtr]), "%$" PRIu64, total);

                        if (strcmp(printNBuffers[currPtr], printNBuffers[otherPtr]))
                        {
                            PrintStatus("%s reads processed", printNBuffers[currPtr]);
                        }

                        printNBufferPtr = otherPtr;
                    }
                }
            }
        } while (readBuffer->size);

        GetNextBuffer_Write(writePool);
        GetNextTransferBuffer(transferBufferPool);
        GetNextTransferBuffer(transferBufferPool);

        WavlTreeFreeze_LowToHigh(barcodeTree);
        
        header = (char *)"Barcode\tCorrect Reads\tCorrected Reads\n";
        if (WriteToLogFile(clearBCLog, header, strlen(header)))
        {
            logError = 1;
            goto End;
        }
        header = (char *)"Barcode\tReads\n";
        if (WriteToLogFile(unclearBCLog, header, strlen(header)))
        {
            logError = 1;
            goto End;
        }
        
        TraverseLinkedList(WavlTreeGetBottom(barcodeTree)->next, wavl_node)
        {
            barcode *barcode = GetBarCodeFromHashTable(barcodeHashTable, &workingSet, node->value - 1);
            u08 lineBuffer[128];

            u08 a = (u08)((barcode->barcode >> 24) & ((1 << 8) - 1));
            u08 c = (u08)((barcode->barcode >> 16) & ((1 << 8) - 1));
            u08 b = (u08)((barcode->barcode >> 8) & ((1 << 8) - 1));
            u08 d = (u08)(barcode->barcode & ((1 << 8) - 1));

            if (barcode->unclear) stbsp_snprintf((char *)lineBuffer, sizeof(lineBuffer), "A%02uC%02uB%02uD%02u\t%u\n", a, c, b, d, barcode->unclear);
            else stbsp_snprintf((char *)lineBuffer, sizeof(lineBuffer), "A%02uC%02uB%02uD%02u\t%u\t%u\n", a, c, b, d, barcode->correct, barcode->corrected);

            if (WriteToLogFile(barcode->unclear ? unclearBCLog : clearBCLog, (char *)lineBuffer, strlen((char *)lineBuffer)))
            {
                logError = 1;
                goto End;
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
