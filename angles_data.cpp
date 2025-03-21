#include "myutils.h"
#include "abcxyz.h"

#define N	790
static double s_AnglesData[N][2] = {
 { 83.4, 159.1}, { 70.3,  37.8}, { 59.7,  30.3}, { 53.2,   9.0}, { 60.6,  29.2},
 { 80.1,   7.9}, { 45.2, 159.2}, {100.0,  84.5}, { 84.4, 124.9}, { 39.8, 141.2},
 { 44.3,  17.1}, { 41.5,  20.0}, { 55.5,  26.0}, { 67.9,  57.8}, { 85.8,   9.7},
 { 42.4,  91.0}, { 85.8, 148.3}, { 70.2,  42.2}, { 66.0,  43.8}, { 69.0,  49.2},
 { 78.2,  14.4}, { 46.1, 162.8}, { 66.7,  17.7}, { 93.0, 160.4}, { 84.6, 128.2},
 { 57.8,  53.3}, { 74.7,  51.5}, { 79.5,   6.6}, { 82.2,  64.2}, { 83.9, 142.9},
 { 53.0,  73.9}, { 57.3,  61.7}, { 86.9,  77.0}, { 88.4, 133.4}, { 93.3, 124.4},
 { 89.4, 131.3}, { 89.3, 125.9}, { 88.9, 131.6}, { 91.8, 125.5}, { 89.8, 130.4},
 { 87.0, 135.3}, { 84.6,  71.9}, { 84.4,   2.9}, { 88.0, 119.9}, { 90.6, 116.7},
 { 60.9,  80.2}, { 67.2,  51.6}, { 92.0, 167.0}, { 86.4, 120.8}, { 55.8,  52.8},
 { 73.9,  42.7}, { 61.4,  29.8}, { 56.1,   2.2}, { 74.7,  20.8}, { 39.5,   2.9},
 { 67.5,  41.9}, { 70.8,  87.8}, { 41.9,   9.5}, { 81.5,  29.0}, { 45.9, 154.2},
 { 62.0,  17.2}, { 61.0,  14.3}, { 55.4,  28.7}, { 40.8,   2.6}, { 37.6, 126.8},
 { 83.3, 122.1}, { 82.0, 126.6}, { 83.1, 138.3}, { 90.5, 130.1}, { 55.7,  38.6},
 { 82.3,  30.0}, { 30.5, 173.3}, { 44.1,  31.1}, { 60.6,  57.4}, { 57.4,  50.9},
 { 55.8,  18.5}, { 91.1,  13.6}, { 53.4, 104.3}, { 40.5,  43.3}, { 80.5, 164.0},
 { 86.4, 129.5}, { 87.0, 141.3}, { 90.5, 119.4}, { 93.5, 122.0}, { 89.4, 129.9},
 { 90.0, 129.7}, { 89.2, 126.4}, { 89.3, 127.8}, { 88.0, 130.5}, { 83.5,  41.1},
 { 77.3,  96.6}, { 37.8,  23.1}, { 85.0, 125.6}, { 62.0,  27.5}, { 49.2,   5.8},
 { 53.4,  13.7}, { 88.3,   3.6}, { 55.4, 130.3}, { 92.7,  36.6}, { 90.0, 117.6},
 { 90.0, 126.2}, { 88.7, 133.5}, { 87.2, 131.3}, { 87.7, 132.8}, { 88.4, 129.1},
 { 84.6, 133.8}, { 90.3, 134.6}, { 85.4, 129.0}, { 76.9, 153.1}, { 90.3,  28.6},
 { 49.2,   5.8}, { 53.4,  13.7}, { 88.3,   3.6}, { 55.4, 130.3}, { 92.7,  36.6},
 { 90.0, 117.6}, { 90.0, 126.2}, { 88.7, 133.5}, { 87.2, 131.3}, { 87.7, 132.8},
 { 88.4, 129.1}, { 84.6, 133.8}, { 90.3, 134.6}, { 85.4, 129.0}, { 76.9, 153.1},
 { 90.3,  28.6}, { 57.8,  37.8}, { 96.0,  50.1}, { 43.8, 177.0}, { 55.4,  27.4},
 { 41.2,   2.9}, { 53.7,  15.3}, { 64.2,  73.6}, { 67.4,  28.8}, { 80.2,   5.2},
 { 38.5, 123.8}, { 50.9,  56.4}, { 68.7,  68.8}, { 80.3,  64.2}, { 42.7,   9.8},
 { 57.3,  12.7}, { 63.3,  21.0}, { 71.5,   3.7}, { 62.9,   7.7}, { 86.5,  24.1},
 { 88.8, 113.6}, { 84.1, 138.5}, { 73.8,  40.5}, { 88.0,  18.7}, { 56.1, 145.3},
 { 37.7,  37.5}, { 85.5,  64.2}, { 52.2, 106.8}, { 63.5,  37.2}, { 42.7,  71.2},
 { 63.1, 120.2}, { 61.8,  64.1}, { 91.7, 125.1}, { 56.1, 145.3}, { 37.7,  37.5},
 { 85.5,  64.2}, { 52.2, 106.8}, { 63.5,  37.2}, { 42.7,  71.2}, { 63.1, 120.2},
 { 61.8,  64.1}, { 91.7, 125.1}, { 49.6, 105.5}, { 41.5,  29.4}, { 45.2,   3.8},
 { 89.6,  94.1}, { 81.0, 121.0}, { 88.1, 145.4}, { 93.5, 121.4}, { 91.4, 132.8},
 { 91.4, 126.3}, { 90.6, 130.8}, { 88.7, 124.7}, { 89.6, 131.2}, { 92.6, 128.7},
 { 87.7, 131.9}, { 88.0, 123.7}, { 87.9, 139.1}, { 89.7, 123.1}, { 89.7, 133.1},
 { 90.3, 127.3}, { 63.1, 120.2}, { 61.8,  64.1}, { 91.7, 125.1}, { 49.6, 105.5},
 { 41.5,  29.4}, { 45.2,   3.8}, { 89.6,  94.1}, { 81.0, 121.0}, { 88.1, 145.4},
 { 93.5, 121.4}, { 91.4, 132.8}, { 91.4, 126.3}, { 90.6, 130.8}, { 75.8, 129.6},
 { 64.1,  35.6}, { 60.5,  89.0}, { 93.9,  46.8}, { 90.4, 119.7}, { 87.9, 127.8},
 { 89.5, 134.7}, { 90.1, 122.5}, { 87.9, 134.1}, { 87.0, 132.1}, { 91.7, 128.7},
 { 77.2, 138.9}, { 72.7, 165.3}, { 74.8, 144.3}, { 88.4,  35.3}, { 91.7, 120.7},
 { 88.1, 104.9}, { 89.1, 123.8}, { 80.6, 131.1}, { 70.8,  48.8}, { 85.3,  43.6},
 { 38.6, 154.6}, { 51.9, 141.3}, { 86.3, 156.9}, { 90.3, 115.5}, { 87.9, 128.2},
 { 93.3, 132.6}, { 87.2, 127.8}, { 96.2, 106.5}, { 85.5,  79.9}, { 60.2,  90.6},
 { 59.5,  80.7}, { 69.1,  93.3}, { 81.9, 140.3}, { 81.5,  76.8}, { 45.2,   9.8},
 { 83.9, 115.1}, { 86.9, 126.9}, { 86.1, 137.5}, { 88.4, 131.0}, { 89.0, 125.6},
 { 89.6, 137.6}, { 86.9, 127.7}, { 85.8, 136.0}, { 90.5, 129.4}, { 85.8, 136.0},
 { 90.5, 129.4}, { 88.1, 129.4}, { 91.6, 128.1}, { 86.5, 127.0}, { 93.4, 123.2},
 { 80.2, 136.5}, { 89.5, 130.1}, { 91.4, 118.7}, { 83.4, 131.1}, { 89.6,  90.8},
 { 48.1, 132.8}, { 83.0,  61.4}, { 86.6,  13.5}, { 91.8, 110.4}, { 84.4, 126.1},
 { 82.4, 130.1}, { 83.4, 135.7}, { 70.4,  80.2}, { 89.1,  11.3}, { 82.6, 134.4},
 { 59.7,  36.8}, { 60.7,  30.1}, { 47.7,  16.9}, { 69.8,  96.4}, { 69.2,  69.0},
 { 67.1,  65.4}, { 61.1,  54.1}, { 62.0,  44.2}, { 73.3, 129.0}, { 42.8,  14.6},
 { 61.9, 102.4}, { 57.3,  17.6}, { 55.1,   6.6}, { 70.7,   3.3}, { 50.3,   3.4},
 { 64.4,  61.8}, { 76.9,   0.4}, { 42.3, 164.9}, { 42.5,   3.4}, { 47.0,   0.2},
 { 65.5,  36.9}, { 74.0,  55.4}, { 85.6,  51.0}, { 82.2, 138.3}, { 71.1, 130.5},
 { 84.3,  20.2}, { 92.4,  92.9}, { 46.2, 134.0}, { 78.4,  22.2}, { 68.4,  28.3},
 { 28.8,  16.2}, { 50.6,  45.4}, { 35.2,  22.0}, { 62.0,  44.2}, { 73.3, 129.0},
 { 42.8,  14.6}, { 61.9, 102.4}, { 57.3,  17.6}, { 55.1,   6.6}, { 70.7,   3.3},
 { 50.3,   3.4}, { 64.4,  61.8}, { 76.9,   0.4}, { 42.3, 164.9}, { 42.5,   3.4},
 { 47.0,   0.2}, { 65.5,  36.9}, { 74.0,  55.4}, { 85.6,  51.0}, { 82.2, 138.3},
 { 71.1, 130.5}, { 84.3,  20.2}, { 92.4,  92.9}, { 46.2, 134.0}, { 78.4,  22.2},
 { 68.4,  28.3}, { 28.8,  16.2}, { 49.5,  65.7}, { 65.2,  27.4}, { 52.2,   2.0},
 { 50.7,   6.7}, { 62.6,  73.3}, { 69.7,   7.6}, { 69.3, 156.2}, { 94.5,  25.0},
 { 70.5, 130.1}, { 90.8, 155.6}, { 80.3,  51.1}, { 69.4,  11.7}, { 43.7,  13.8},
 { 61.4,  49.2}, { 61.3,   6.9}, { 58.8,   0.3}, { 63.0,  50.6}, { 38.2,  18.2},
 { 51.0,  18.2}, { 34.3,  82.1}, { 36.4,  82.0}, { 58.8,  68.2}, { 66.7,  37.2},
 { 88.6,  15.4}, { 79.0, 126.2}, { 92.5, 152.4}, { 74.6,  44.3}, { 49.9,  17.5},
 { 83.4,   4.9}, { 45.7, 141.0}, { 92.1,  62.3}, { 89.2, 121.9}, { 91.3, 133.9},
 { 86.4, 127.7}, { 87.6, 135.6}, { 91.5, 126.2}, { 90.9, 112.0}, { 81.0, 108.1},
 { 89.3, 145.0}, { 93.9, 122.0}, { 52.2,  47.8}, { 89.8, 113.2}, { 90.6,  41.5},
 { 86.6, 122.4}, { 90.8, 131.4}, { 85.0, 126.8}, { 85.0,  76.9}, { 84.7,  58.3},
 { 59.2, 123.0}, { 44.8,  17.5}, { 72.4,  42.8}, { 47.3,   1.7}, { 73.6,  24.6},
 { 43.7, 119.2}, { 94.7,  49.9}, { 89.5, 120.2}, { 86.7, 134.3}, { 92.2, 125.8},
 { 87.3, 135.4}, { 88.9, 128.8}, { 88.6, 133.7}, { 90.3, 132.1}, { 90.7, 131.3},
 { 85.0, 127.5}, { 66.8, 152.9}, { 70.5,  76.6}, { 75.6, 179.6}, { 91.2,  43.8},
 { 85.2, 119.2}, { 91.5, 125.7}, { 91.0, 131.6}, { 89.6, 129.8}, { 85.6, 137.0},
 { 88.2, 131.2}, { 90.0, 136.8}, { 91.2, 115.2}, { 43.2,  45.4}, { 77.9,  77.5},
 { 82.5,  95.9}, { 44.6,  66.9}, { 87.4,  35.2}, { 69.9,   5.0}, { 78.6,   2.8},
 { 41.5, 147.7}, { 46.7,   3.0}, { 50.2,   3.7}, { 44.7,  21.2}, { 63.6,  10.2},
 { 18.7,   5.1}, { 72.7, 147.9}, { 79.6, 170.7}, { 66.1, 127.4}, { 40.7, 155.2},
 { 81.4,   8.1}, { 84.7,  96.6}, { 63.3,  13.5}, { 38.9, 166.2}, { 66.7,   9.1},
 { 44.2,   3.4}, { 56.5,  24.2}, { 43.4,   3.0}, { 86.0, 167.4}, { 60.7,  28.8},
 { 62.3,  68.8}, { 85.2,  62.9}, { 89.3, 127.1}, { 87.7, 129.5}, { 89.4, 130.8},
 { 89.5, 127.4}, { 88.6, 135.6}, { 88.7, 126.0}, { 88.8, 126.4}, { 79.6, 138.7},
 { 79.7, 144.8}, { 55.8, 175.9}, { 84.0, 125.0}, { 84.3, 132.6}, { 86.4, 138.6},
 { 85.3, 125.7}, { 78.2,  68.7}, { 66.3,  20.8}, { 79.0,  34.1}, { 27.7, 140.5},
 { 49.5,  77.5}, { 73.8,  34.3}, { 78.2,  68.7}, { 66.3,  20.8}, { 79.0,  34.1},
 { 27.7, 140.5}, { 49.5,  77.5}, { 73.8,  34.3}, { 50.3,  22.5}, { 63.5,  77.5},
 { 72.9, 116.8}, { 84.2,  18.9}, { 54.6,  82.9}, { 64.0,  95.2}, { 41.9,   8.4},
 { 62.4,  88.2}, { 68.0,  62.2}, { 66.0,  24.2}, { 70.1,  51.4}, { 83.1,  57.8},
 { 57.5,  61.7}, { 37.8,  47.8}, { 94.3, 100.4}, { 89.5, 117.7}, { 33.4,  97.2},
 { 47.9, 109.7}, { 85.5,  94.3}, { 88.9, 127.4}, { 89.0, 133.7}, { 87.4, 128.9},
 { 89.7, 132.9}, { 90.3, 127.0}, { 84.9, 131.1}, { 88.7, 132.6}, { 85.7, 132.5},
 { 84.7, 134.9}, { 85.7, 133.1}, { 91.0, 132.2}, { 89.6, 117.4}, { 60.4,  27.3},
 { 77.0,  17.4}, { 84.9,   4.7}, { 89.0, 111.4}, { 87.0, 159.9}, { 67.4,  62.1},
 { 50.2,  29.0}, { 44.6,   3.5}, { 91.5, 101.9}, { 89.6, 112.5}, { 82.4, 119.2},
 { 45.6, 102.0}, { 76.5,  57.0}, { 87.5,  59.5}, { 93.5,  19.5}, { 90.5, 103.8},
 { 68.6, 110.9}, { 88.2,  15.9}, { 86.9, 118.5}, { 85.9, 136.5}, { 91.2, 122.6},
 { 75.8, 147.2}, { 77.0, 143.7}, { 48.6, 110.5}, { 60.7, 101.6}, { 74.5,  40.2},
 { 90.2, 157.5}, { 83.5, 124.5}, { 64.0,  33.1}, { 75.0,  77.4}, { 65.0,  43.4},
 { 61.1,  47.8}, { 85.0,  64.1}, { 88.6,   3.1}, { 89.7, 109.8}, { 61.4,  88.7},
 { 89.4, 115.3}, { 76.8, 132.1}, { 40.1, 107.7}, { 57.0, 153.7}, { 49.3,  35.6},
 { 71.0,  75.4}, { 78.3,  35.6}, { 89.5,   6.3}, { 92.9, 108.8}, { 79.1, 111.7},
 { 85.7, 147.7}, { 86.9, 123.4}, { 85.0, 143.0}, { 88.7, 128.7}, { 88.7, 128.4},
 { 88.8, 131.7}, { 84.1, 132.6}, { 89.0, 133.7}, { 91.2, 130.4}, { 86.7, 120.4},
 { 90.6,  64.0}, { 82.4, 116.0}, { 68.2,  51.2}, { 50.3,  37.6}, { 55.9, 120.4},
 { 73.1, 112.3}, { 32.1,  16.8}, { 74.1,  47.3}, { 88.0, 149.1}, { 84.3, 137.7},
 { 70.5,  62.9}, { 65.8,  62.1}, { 57.6,  58.7}, { 60.9,  70.3}, { 73.1,  71.1},
 { 91.7,  53.4}, { 83.2,  37.4}, { 48.2, 137.4}, { 45.8,   1.2}, { 75.5,   5.1},
 { 61.4,   8.5}, { 64.7,  74.3}, { 67.2, 115.7}, { 88.6, 171.6}, { 33.5,  20.5},
 { 81.1, 128.7}, { 91.2,  12.6}, { 87.5,  30.3}, { 71.4,  28.6}, { 87.1,  44.4},
 { 48.2, 137.4}, { 45.8,   1.2}, { 75.5,   5.1}, { 61.4,   8.5}, { 64.7,  74.3},
 { 67.2, 115.7}, { 88.6, 171.6}, { 33.5,  20.5}, { 81.1, 128.7}, { 91.2,  12.6},
 { 87.5,  30.3}, { 71.4,  28.6}, { 87.1,  44.4}, { 90.7, 115.9}, { 77.0, 139.2},
 { 71.1,  39.0}, { 56.4,  83.7}, { 54.9,  51.1}, { 40.5,  10.8}, { 77.8,  64.5},
 { 48.3,  45.2}, { 89.8,  82.3}, { 61.0,  22.3}, { 48.3,  16.9}, { 74.3,  41.1},
 { 87.9,   3.6}, { 90.0,  33.0}, { 87.6, 127.1}, { 56.4,  47.9}, { 86.3,  23.1},
 { 74.4,  30.7}, { 86.5,  30.2}, { 79.1, 129.3}, { 73.1,  56.2}, { 74.8,  33.2},
 { 69.8,  42.3}, { 69.4,  58.0}, { 86.5,  30.2}, { 79.1, 129.3}, { 73.1,  56.2},
 { 74.8,  33.2}, { 69.8,  42.3}, { 69.4,  58.0}, { 72.1,  76.5}, { 52.5,  17.8},
 { 37.5,  17.6}, { 74.5,  95.2}, { 35.0,  87.1}, { 81.6, 154.1}, { 45.0, 140.3},
 { 90.2,  52.9}, { 43.9, 146.6}, { 92.0,  57.2}, { 40.5, 154.6}, { 58.8, 128.6},
 { 68.6,  67.6}, { 67.9,  13.7}, { 73.0,  13.5}, { 68.4,  10.9}, { 46.4,   5.2},
 { 94.4, 143.5}, { 65.5,   1.0}, { 50.4, 154.7}, { 66.3, 133.1}, { 81.4,  57.4},
 { 76.9,  63.4}, { 59.5,  45.1}, { 86.9,  71.6}, { 87.5, 129.9}, { 88.8, 131.3},
 { 89.5, 129.5}, { 89.7, 135.0}, { 85.4, 124.4}, { 83.8, 139.5}, { 91.1, 130.9},
 { 87.8, 122.3}, { 93.5, 117.1}, { 90.7, 130.5}, { 88.3, 128.3}, { 86.3, 134.6},
 { 90.7, 130.5}, { 88.3, 128.3}, { 86.3, 134.6}, { 90.4, 128.8}, { 86.2, 132.9},
 { 56.7, 152.2}, { 88.2, 176.0}, { 85.1,  33.7}, { 72.5,  11.6}, { 57.1,  53.9},
 { 79.9,  72.2}, { 40.6, 105.7}, { 60.2, 172.8}, { 86.2, 132.9}, { 56.7, 152.2},
 { 88.2, 176.0}, { 85.1,  33.7}, { 72.5,  11.6}, { 57.1,  53.9}, { 79.9,  72.2},
 { 40.6, 105.7}, { 60.2, 172.8}, { 85.3, 153.5}, { 92.6, 122.2}, { 87.8, 129.5},
 { 90.2, 130.7}, { 89.5, 130.2}, { 89.7, 126.9}, { 87.3, 136.7}, { 67.2,  40.6},
 { 85.8, 174.4}, { 86.6,  38.3}, { 85.4, 123.8}, { 77.8, 147.8}, { 71.3,  75.1},
 { 66.7,  52.8}, { 55.7, 176.1}, { 51.1, 165.4}, { 51.5, 175.4}, { 40.8,  96.5},
 { 83.6,  79.0}, { 44.3, 168.4}, { 50.6,  32.0}, { 64.0,  37.6}, { 86.9,  27.1},
 { 46.0, 160.9}, { 67.3,  19.7}, { 55.7,  10.8}, { 71.7,  42.4}, { 61.6,   9.3},
 { 64.5,  34.0}, { 61.1,   1.9}, { 87.9, 156.1}, { 84.8, 132.7}, { 60.1,  69.3},
 { 67.7,  41.9}, { 55.4,   7.1}, { 61.9,  49.7}, { 62.8,  20.8}, { 66.5,  31.6},
 { 52.8,  17.3}, { 60.0,  15.2}, { 60.8,  47.5}, { 65.4,  62.2}, { 58.4,  26.9},
 { 53.6,   3.5}, { 69.5,  43.3}, { 63.1,   1.7}, { 41.9,  16.7}, { 91.8,  39.0},
 { 53.6,  80.7}, { 59.6,  85.1}, { 88.7,  56.5}, { 95.7, 109.2}, { 73.3, 155.0},
 { 80.1, 136.9}, { 85.9, 121.5}, { 88.8, 135.5}, { 88.9, 129.3}, { 89.9, 131.0},
 { 88.0, 126.1}, { 86.1, 141.4}, { 75.5, 135.4}, { 76.3, 115.9}, { 86.8,  21.9},
 { 81.7, 128.4}, { 86.3, 127.4}, { 60.0,  89.5}, { 66.7,  57.1}, { 75.7,  40.6},
 { 54.0,  15.5}, { 84.7, 163.5}, { 84.0, 131.9}, { 46.0,  86.2}, { 74.1,  56.2},
 { 54.2,  19.2}, { 83.1,  44.9}, { 79.7,  21.2}, { 90.5,  12.0}, { 91.0, 106.3},
 { 88.9, 120.8}, { 83.9, 100.5}, { 84.2,  49.4}, { 49.5,  60.1}, { 78.2,  85.8},
 { 91.8,  36.5}, { 97.0, 114.0}, { 82.3, 136.2}, { 46.9,  87.7}, { 48.3,   0.3},
 { 54.2,  25.4}, { 72.2,  14.7}, { 54.0,  15.5}, { 84.7, 163.5}, { 84.0, 131.9},
 { 46.0,  86.2}, { 74.1,  56.2}, { 54.2,  19.2}, { 83.1,  44.9}, { 79.7,  21.2},
 { 90.5,  12.0}, { 91.0, 106.3}, { 88.9, 120.8}, { 83.9, 100.5}, { 84.2,  49.4},
 { 49.5,  60.1}, { 90.7,   9.3}, { 89.1, 112.8}, { 62.7,  54.3}, { 79.0,   8.5},
 { 83.7,  63.0}, { 88.6, 133.4}, { 86.3, 133.1}, { 88.3, 131.0}, { 83.2, 119.9},
 { 82.5,  98.7}, { 91.9, 126.8}, { 86.7, 129.2}, { 92.7, 137.5}, { 40.6,  74.5},
};

void GetRandomAnglePair_Degrees(double &theta_deg, double &phi_deg)
	{
	uint i = randu32()%N;
	theta_deg = s_AnglesData[i][0];
	phi_deg = s_AnglesData[i][1];
	}

void GetRandomAnglePair_Radians(double &theta_rad, double &phi_rad)
	{
	double theta_deg, phi_deg;
	GetRandomAnglePair_Degrees(theta_deg, phi_deg);
	theta_rad = radians(theta_deg);
	phi_rad = radians(phi_deg);
	}
