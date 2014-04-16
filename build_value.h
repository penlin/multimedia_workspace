#ifndef __BUILD_VALUE_H
#define __BUILD_VALUE_H

/**#_______  _______  _        _______ _________ _______  _       _________              _______  _                 _______
#  (  ____ \(  ___  )( (    /|(  ____ \\__   __/(  ___  )( (    /|\__   __/    |\     /|(  ___  )( \      |\     /|(  ____ \
#  | (    \/| (   ) ||  \  ( || (    \/   ) (   | (   ) ||  \  ( |   ) (       | )   ( || (   ) || (      | )   ( || (    \/
#  | |      | |   | ||   \ | || (_____    | |   | (___) ||   \ | |   | |       | |   | || (___) || |      | |   | || (__
#  | |      | |   | || (\ \) |(_____  )   | |   |  ___  || (\ \) |   | |       ( (   ) )|  ___  || |      | |   | ||  __)
#  | |      | |   | || | \   |      ) |   | |   | (   ) || | \   |   | |        \ \_/ / | (   ) || |      | |   | || (
#  | (____/\| (___) || )  \  |/\____) |   | |   | )   ( || )  \  |   | |         \   /  | )   ( || (____/\| (___) || (____/\
#  (_______/(_______)|/    )_)\_______)   )_(   |/     \||/    )_)   )_(          \_/   |/     \|(_______/(_______)(_______/
# **/
#define PXL 8

// decode algorithm function name
#define __DIRECT__      direct_system
#define __INTRA__       intra_system
#define __INTER__       inter_system
#define __INTER_INTRA__ inter_intra

// output sequence type
#define __YUV420__    1
#define __Y__         2

// BCJR decoding algorithm
#define __MAP__     1
#define __LOG_MAP__ 2


/**#______   _______  _______ _________ _______         _________ _        _______  _______  _______  _______  _______ __________________ _______  _
#  (  ___ \ (  ___  )(  ____ \\__   __/(  ____ \        \__   __/( (    /|(  ____ \(  ___  )(  ____ )(       )(  ___  )\__   __/\__   __/(  ___  )( (    /|
#  | (   ) )| (   ) || (    \/   ) (   | (    \/           ) (   |  \  ( || (    \/| (   ) || (    )|| () () || (   ) |   ) (      ) (   | (   ) ||  \  ( |
#  | (__/ / | (___) || (_____    | |   | |                 | |   |   \ | || (__    | |   | || (____)|| || || || (___) |   | |      | |   | |   | ||   \ | |
#  |  __ (  |  ___  |(_____  )   | |   | |                 | |   | (\ \) ||  __)   | |   | ||     __)| |(_)| ||  ___  |   | |      | |   | |   | || (\ \) |
#  | (  \ \ | (   ) |      ) |   | |   | |                 | |   | | \   || (      | |   | || (\ (   | |   | || (   ) |   | |      | |   | |   | || | \   |
#  | )___) )| )   ( |/\____) |___) (___| (____/\        ___) (___| )  \  || )      | (___) || ) \ \__| )   ( || )   ( |   | |   ___) (___| (___) || )  \  |
#  |/ \___/ |/     \|\_______)\_______/(_______/        \_______/|/    )_)|/       (_______)|/   \__/|/     \||/     \|   )_(   \_______/(_______)|/    )_)
# **/

// basic decoding information
#define Niter       3                  // #iteration
#define __HEIGHT    288
#define __WIDTH     352
#define __SKIP      0
#define __FRAME     300
#define __SNR       4
#define __SNR_S     0
#define __SNR_E     4


#define __SEQ_DIR   "sequence/"
#define __FOREMAN   __SEQ_DIR  "foreman_cif.yuv"
#define __STEFAN    __SEQ_DIR  "stefan_cif.yuv"
#define __HALL      __SEQ_DIR  "hall_cif.yuv"
#define __AKIYO     __SEQ_DIR  "akiyo_cif.yuv"
#define __TAG__     "foreman"

/**#______   _______  _______  _______  ______   _______        _______  _______  _       _________ _______  _______  _
#  (  __  \ (  ____ \(  ____ \(  ___  )(  __  \ (  ____ \      (  ____ \(  ___  )( (    /|\__   __/(  ____ )(  ___  )( \
#  | (  \  )| (    \/| (    \/| (   ) || (  \  )| (    \/      | (    \/| (   ) ||  \  ( |   ) (   | (    )|| (   ) || (
#  | |   ) || (__    | |      | |   | || |   ) || (__          | |      | |   | ||   \ | |   | |   | (____)|| |   | || |
#  | |   | ||  __)   | |      | |   | || |   | ||  __)         | |      | |   | || (\ \) |   | |   |     __)| |   | || |
#  | |   ) || (      | |      | |   | || |   ) || (            | |      | |   | || | \   |   | |   | (\ (   | |   | || |
#  | (__/  )| (____/\| (____/\| (___) || (__/  )| (____/\      | (____/\| (___) || )  \  |   | |   | ) \ \__| (___) || (____/\
#  (______/ (_______/(_______/(_______)(______/ (_______/      (_______/(_______)|/    )_)   )_(   |/   \__/(_______)(_______/
# **/

// control the decode algorithm
#define __ALGO__    __INTRA__
#define __SEQ__     __HALL

#define __RANDOM__      1               // total control for the interleave and AWGN
#define __INTERLEAVE__  1&&__RANDOM__   // generate random order map
#define __NOISE__       1&&__RANDOM__   // generate gaussian noise

#define __DEBUG__       0              // debug msg
#define __STATUS__      0&&__DEBUG__    // msg for status right now
#define __PROGRESS__    1&&__DEBUG__    // msg for the progress
#define __BETA__        0&&__DEBUG__    // msg for beta msg
#define __PSNR__        1&&__DEBUG__    // msg for cPSNR in procedure

#define __OUTPUT_SEQ__      0           // control if output the decoded sequence
#define __OUTPUT_TYPE__     __Y__*__OUTPUT_SEQ__

#define __EXIT_INFO__       0

#define __BCJR__    __MAP__

#endif // __BUILD_VALUE_H
