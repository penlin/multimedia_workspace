function p = ps_bp6()
Ia=[0.413264,0.412560,0.413665,0.416202,0.414647,0.416909,0.419268,0.419555,0.417903,0.419516,0.421047,0.420263,0.422761,0.421889,0.424364,0.424596,0.423152,0.425218,0.427780,0.424386,0.426287,0.430860,0.429073,0.429086,0.432267,0.431175,0.431657,0.434570,0.432995,0.435342,0.435802,0.435036,0.434943,0.437517,0.434790,0.438367,0.438409,0.437502,0.441413,0.444307,0.442325,0.443950,0.444639,0.443934,0.446746,0.447743,0.449671,0.447856,0.447212,0.449925,0.452177,0.451883,0.449287,0.453344,0.454169,0.454067,0.454401,0.454514,0.454322,0.458290,0.455367,0.458570,0.459372,0.458512,0.461129,0.462295,0.461125,0.464582,0.462351,0.462203,0.465426,0.465192,0.467340,0.466758,0.468567,0.470095,0.468535,0.470938,0.472815,0.473870,0.471629,0.472675,0.474937,0.478789,0.477361,0.475717,0.476410,0.478191,0.479021,0.481728,0.480286,0.481424,0.482107,0.482709,0.482373,0.486409,0.486117,0.486850,0.485443,0.486664,0.489714,0.489682,0.489771,0.487775,0.490947,0.491290,0.492126,0.493066,0.495014,0.497713,0.498007,0.498142,0.495437,0.497808,0.501293,0.502193,0.501684,0.501604,0.502828,0.501675,0.504943,0.502720,0.507075,0.505446,0.507611,0.510201,0.507742,0.509470,0.512262,0.511060,0.515047,0.514962,0.512844,0.514741,0.516787,0.518224,0.514795,0.516097,0.516388,0.517186,0.522358,0.520100,0.523800,0.523271,0.524247,0.524526,0.526889,0.525306,0.527477,0.532101,0.524591,0.529908,0.526816,0.531374,0.532381,0.534741,0.533977,0.535589,0.533857,0.534411,0.534295,0.539305,0.539014,0.538093,0.541355,0.541981,0.538554,0.541142,0.540432,0.542346,0.546465,0.545441,0.546462,0.547716,0.548888,0.547720,0.548008,0.548430,0.551028,0.551280,0.555450,0.553231,0.552924,0.558657,0.556853,0.555718,0.558695,0.557821,0.559819,0.558386,0.562848,0.562238,0.562961,0.562625,0.564042,0.566171,0.566159,0.566181,0.568634,0.566971,0.571123,0.570631,0.569897,0.572901,0.570925,0.572032,0.571640,0.576861,0.576529,0.576409,0.579021,0.577507,0.575882,0.579536,0.581420,0.580310,0.582217,0.584480,0.584399,0.583067,0.590164,0.584726,0.588598,0.587389,0.587850,0.588531,0.589774,0.593462,0.592174,0.594337,0.595074,0.594361,0.594454,0.595133,0.595710,0.596539,0.600346,0.599824,0.601556,0.601826,0.599432,0.603688,0.604123,0.605089,0.606776,0.605767,0.607780,0.606402,0.608436,0.611367,0.610364,0.609674,0.610858,0.614296,0.614908,0.614552,0.615816,0.614893,0.616241,0.616708,0.617821,0.618572,0.619552,0.622647,0.621920,0.621358,0.623250,0.626133,0.625787,0.627477,0.627740,0.628132,0.627974,0.630884,0.631451,0.628867,0.631765,0.633223,0.630649,0.636168,0.636229,0.637943,0.637245,0.637069,0.635433,0.640761,0.642283,0.640922,0.640315,0.643059,0.643681,0.643615,0.642906,0.645544,0.648381,0.647471,0.649317,0.649259,0.652598,0.650860,0.651593,0.652502,0.652855,0.654880,0.654225,0.656100,0.655801,0.661490,0.660245,0.655172,0.661422,0.659042,0.663684,0.660430,0.662959,0.663617,0.665141,0.664075,0.664761,0.668508,0.668646,0.669318,0.672973,0.669749,0.672699,0.671875,0.673362,0.674174,0.676677,0.676306,0.678395,0.675660,0.677907,0.680906,0.679464,0.681702,0.677856,0.682287,0.685136,0.682765,0.684347,0.685337,0.687323,0.686359,0.688840,0.690567,0.691090,0.692245,0.690877,0.693186,0.693869,0.691721,0.694991,0.697290,0.699603,0.697092,0.699193,0.699966,0.700354,0.698642,0.701535,0.702931,0.701637,0.703755,0.704048,0.708932,0.706483,0.708787,0.708575,0.707442,0.708434,0.710284,0.710532,0.711338,0.712805,0.716549,0.714365,0.714197,0.715931,0.716815,0.716549,0.719072,0.719748,0.721060,0.718258,0.722564,0.724142,0.724615,0.724201,0.725486,0.725861,0.726653,0.728349,0.730013,0.728127,0.729585,0.730077,0.731809,0.731716,0.732023,];
Ie=[0.025475,0.025501,0.025686,0.025744,0.024791,0.026416,0.026447,0.027736,0.025879,0.026483,0.027056,0.026366,0.026026,0.026677,0.026555,0.026128,0.025848,0.026480,0.027937,0.026051,0.027275,0.026803,0.027298,0.026243,0.027279,0.028076,0.026610,0.028368,0.027568,0.027105,0.028363,0.027167,0.027443,0.028205,0.028208,0.028506,0.028372,0.029196,0.029279,0.027543,0.027670,0.029198,0.028342,0.029162,0.028494,0.028424,0.029790,0.028353,0.028904,0.029048,0.030266,0.029252,0.029195,0.028791,0.030589,0.029730,0.029494,0.029785,0.030087,0.030157,0.030451,0.030659,0.029371,0.030527,0.030128,0.032718,0.029824,0.032424,0.030770,0.031563,0.029552,0.030674,0.031140,0.030986,0.032559,0.030881,0.031187,0.031830,0.031221,0.031677,0.032073,0.030955,0.032980,0.033204,0.033722,0.032341,0.032022,0.032232,0.030928,0.031878,0.031580,0.033092,0.033036,0.032029,0.031950,0.033483,0.033624,0.033663,0.032475,0.033086,0.034372,0.032778,0.034206,0.033363,0.033639,0.034118,0.033986,0.035906,0.034779,0.035798,0.034620,0.034780,0.033789,0.034598,0.034430,0.035658,0.034423,0.034254,0.034284,0.034509,0.035145,0.035222,0.036289,0.035369,0.036164,0.035150,0.035990,0.036395,0.035279,0.036230,0.035181,0.035441,0.034879,0.035931,0.036343,0.036689,0.037169,0.037443,0.037278,0.037390,0.038847,0.036684,0.037955,0.037438,0.036309,0.037327,0.038052,0.036969,0.038067,0.037209,0.038496,0.037038,0.037947,0.037810,0.038806,0.038423,0.039751,0.038972,0.037783,0.038006,0.038268,0.039516,0.038410,0.038697,0.038916,0.039695,0.039482,0.038923,0.038717,0.038726,0.039960,0.040607,0.040146,0.039607,0.040132,0.040615,0.040370,0.040754,0.040140,0.039737,0.042094,0.040623,0.040977,0.040356,0.041222,0.041307,0.040259,0.041133,0.040687,0.040734,0.042361,0.042583,0.041971,0.042667,0.042215,0.042328,0.042492,0.040616,0.042284,0.042447,0.043597,0.042588,0.042142,0.043366,0.041453,0.042173,0.042361,0.043085,0.043014,0.043100,0.043721,0.041985,0.043459,0.043402,0.044287,0.044777,0.042202,0.044220,0.045065,0.043011,0.045855,0.043393,0.044929,0.044965,0.044951,0.045882,0.044497,0.044334,0.044300,0.045656,0.046146,0.044674,0.044436,0.044756,0.045006,0.045484,0.045921,0.045327,0.044622,0.047088,0.046089,0.045977,0.046428,0.046210,0.045823,0.046674,0.046537,0.046873,0.046795,0.046146,0.046896,0.047334,0.047008,0.047615,0.047677,0.047187,0.047713,0.046557,0.048094,0.047681,0.048887,0.047863,0.048413,0.048844,0.047485,0.047889,0.048425,0.049698,0.048951,0.048016,0.049403,0.048463,0.048323,0.049005,0.049577,0.049257,0.049141,0.050244,0.049505,0.049277,0.049954,0.051625,0.049354,0.049937,0.050496,0.050225,0.051582,0.050987,0.049294,0.050111,0.051258,0.050798,0.049508,0.051115,0.051188,0.050689,0.050831,0.051934,0.051331,0.050496,0.051744,0.051709,0.051636,0.051770,0.050767,0.053256,0.052190,0.052611,0.052458,0.052255,0.052029,0.052953,0.052209,0.052058,0.051675,0.051086,0.052870,0.053024,0.053584,0.053921,0.053325,0.052853,0.054995,0.053867,0.053685,0.054584,0.053223,0.054203,0.053827,0.054652,0.054118,0.054045,0.053698,0.053989,0.055007,0.054425,0.053856,0.055835,0.055191,0.054477,0.055805,0.054289,0.055794,0.054599,0.055529,0.056265,0.055605,0.055361,0.056046,0.056122,0.056475,0.055814,0.056069,0.055644,0.056640,0.055928,0.057956,0.057156,0.055759,0.057104,0.056351,0.057488,0.057635,0.057348,0.056584,0.057305,0.058018,0.057244,0.057135,0.056461,0.056412,0.057482,0.056733,0.058058,0.058427,0.057630,0.058284,0.057737,0.059563,0.058170,0.059807,0.058395,0.058062,0.059937,0.058692,0.059131,0.060753,0.059666,0.060208,0.060136,0.059818,0.060104,0.058106,0.060111,0.060399,0.060292,0.059564,0.060249,0.060835,0.060554,];
p=polyfit(Ia,Ie,3);
xi = 0:0.01:1;
yi = polyval(p,xi);
plot(Ie, Ia,'o','LineWidth',2,'Color',[0 0 0]);
hold on;
plot(yi,xi,'LineWidth',1);
axis([0,1,0,1]);
grid on;
hold off;