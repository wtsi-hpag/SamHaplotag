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

#include "Header.h"

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wreserved-id-macro"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wextra-semi-stmt"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wconditional-uninitialized"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#pragma GCC diagnostic ignored "-Wpadded"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#define STB_SPRINTF_IMPLEMENTATION
#include "stb_sprintf.h"
#pragma clang diagnostic pop

#include "WAVLTree.cpp"

global_variable
char
Message_Buffer[1024];

#define PrintError(message, ...) do \
{ \
    stbsp_snprintf(Message_Buffer, 512, message, ##__VA_ARGS__); \
    stbsp_snprintf(Message_Buffer + 512, 512, "[" ProgramName " Error] :: %s\n", Message_Buffer); \
    fprintf(stderr, "%s", Message_Buffer + 512); \
} while (0)

#define PrintStatus(message, ...) do \
{ \
    stbsp_snprintf(Message_Buffer, 512, message, ##__VA_ARGS__); \
    stbsp_snprintf(Message_Buffer + 512, 512, "[" ProgramName " Status] :: %s\n", Message_Buffer); \
    fprintf(stderr, "%s", Message_Buffer + 512); \
} while (0)

#define PrintWarning(message, ...) do \
{ \
    stbsp_snprintf(Message_Buffer, 512, message, ##__VA_ARGS__); \
    stbsp_snprintf(Message_Buffer + 512, 512, "[" ProgramName " Warning] :: %s\n", Message_Buffer); \
    fprintf(stderr, "%s", Message_Buffer + 512); \
} while (0)

struct
buffer
{
    u08 *buffer;
    u64 size;
};

struct
buffer_pool
{
    thread_pool *pool;
    s32 handle;
    u32 bufferPtr;
    buffer *buffers[2];
};

global_function
buffer_pool *
CreatePool(memory_arena *arena)
{
    buffer_pool *pool = PushStructP(arena, buffer_pool);
    pool->pool = ThreadPoolInit(arena, 1);

#define BufferSize MegaByte(16)
    pool->bufferPtr = 0;
    pool->buffers[0] = PushStructP(arena, buffer);
    pool->buffers[0]->buffer = PushArrayP(arena, u08, BufferSize);
    pool->buffers[0]->size = 0;
    pool->buffers[1] = PushStructP(arena, buffer);
    pool->buffers[1]->buffer = PushArrayP(arena, u08, BufferSize);
    pool->buffers[1]->size = 0;

    return(pool);
}

global_function
void
FillBuffer(void *in)
{
    buffer_pool *pool = (buffer_pool *)in;
    buffer *buffer = pool->buffers[pool->bufferPtr];
    buffer->size = (u64)read(pool->handle, buffer->buffer, BufferSize);
}

global_variable
u08
Global_Write_Error = 0;

global_function
void
OutputBuffer(void *in)
{
    buffer_pool *pool = (buffer_pool *)in;
    buffer *buffer = pool->buffers[pool->bufferPtr];
    if ((u64)write(pool->handle, buffer->buffer, buffer->size) != buffer->size) Global_Write_Error = 1;
}

global_function
buffer *
GetNextBuffer_Read(buffer_pool *pool)
{
    FenceIn(ThreadPoolWait(pool->pool));
    buffer *buffer = pool->buffers[pool->bufferPtr];
    pool->bufferPtr = (pool->bufferPtr + 1) & 1;
    ThreadPoolAddTask(pool->pool, FillBuffer, pool);
    return(buffer);
}

global_function
buffer *
GetNextBuffer_Write(buffer_pool *pool)
{
    FenceIn(ThreadPoolWait(pool->pool));
    buffer *buffer = pool->buffers[pool->bufferPtr];
    pool->bufferPtr = (pool->bufferPtr + 1) & 1;
    ThreadPoolAddTask(pool->pool, OutputBuffer, pool);
    buffer->size = 0;
    return(buffer);
}

global_function
u08
WriteToLogFile(s32 handle, void *buffer, u64 size)
{
    return((u64)write(handle, buffer, size) != size);
}

