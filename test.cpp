#include <stdio.h>
#include <stdlib.h>
#include "data_alloc.h"
#include "channel_code_utils2.h"

int main(){
    initMemMgr(65536);

    uint8* x = MALLOC(uint8,10);



    freeMemMgr();
}
