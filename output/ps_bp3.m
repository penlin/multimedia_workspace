function p = ps_bp3()
Ia=[0.412936,0.414355,0.414744,0.415934,0.418605,0.416838,0.418712,0.420021,0.419271,0.418093,0.422983,0.422431,0.424136,0.423092,0.422578,0.424177,0.427400,0.426332,0.426937,0.427193,0.425487,0.429432,0.430720,0.430219,0.432489,0.431148,0.433650,0.436226,0.434060,0.434216,0.436320,0.434385,0.436986,0.436501,0.438718,0.437154,0.440757,0.437853,0.442579,0.441108,0.442749,0.443133,0.444108,0.445642,0.446776,0.444537,0.447407,0.451273,0.446774,0.446949,0.446833,0.449556,0.451743,0.452262,0.453715,0.455685,0.455776,0.456573,0.457556,0.458632,0.459194,0.456217,0.460226,0.459838,0.461354,0.462368,0.462492,0.464647,0.464154,0.464357,0.466862,0.465987,0.466632,0.466737,0.469925,0.466095,0.469870,0.468520,0.470805,0.473412,0.473241,0.474289,0.474614,0.475380,0.474067,0.478896,0.475530,0.481425,0.480231,0.481728,0.480758,0.482384,0.483374,0.484515,0.482777,0.483211,0.486449,0.488691,0.484626,0.486531,0.488298,0.486888,0.492746,0.493455,0.490979,0.490468,0.493405,0.491608,0.495474,0.495134,0.497358,0.498989,0.498923,0.500615,0.499867,0.500830,0.500364,0.502863,0.505101,0.502704,0.504842,0.504169,0.503705,0.503343,0.504940,0.505798,0.511374,0.507723,0.512044,0.512460,0.510610,0.514045,0.514615,0.513568,0.517525,0.515985,0.518216,0.520659,0.519028,0.517976,0.519019,0.520484,0.521863,0.521685,0.524859,0.524113,0.524887,0.525657,0.527034,0.525719,0.528629,0.528295,0.526837,0.530522,0.534391,0.530415,0.532290,0.533361,0.535233,0.537182,0.538793,0.538215,0.537658,0.539666,0.539308,0.536622,0.539863,0.542554,0.542724,0.542615,0.545086,0.545067,0.546062,0.543900,0.546757,0.549538,0.549160,0.550352,0.551049,0.551746,0.552948,0.553681,0.553135,0.557119,0.553566,0.558346,0.557703,0.556253,0.558446,0.560183,0.557446,0.560190,0.561529,0.563835,0.564675,0.566589,0.565103,0.568060,0.566734,0.569790,0.567937,0.569276,0.570722,0.573403,0.572709,0.574315,0.569969,0.573710,0.577161,0.575075,0.573734,0.576417,0.577381,0.577182,0.579926,0.584072,0.582328,0.584608,0.585850,0.582690,0.585866,0.589430,0.587059,0.589273,0.589031,0.588186,0.590985,0.591515,0.592029,0.591644,0.596354,0.595992,0.594599,0.592334,0.597356,0.595162,0.596943,0.599304,0.597316,0.599787,0.600709,0.603198,0.603253,0.603171,0.604862,0.607348,0.604732,0.607917,0.610102,0.612210,0.608814,0.611824,0.613362,0.614996,0.612041,0.614430,0.615823,0.616127,0.616742,0.618760,0.617444,0.619216,0.619193,0.624449,0.622960,0.626163,0.623767,0.626479,0.624940,0.624774,0.626379,0.627287,0.628676,0.627529,0.631297,0.630822,0.630526,0.634697,0.635135,0.633923,0.635989,0.634889,0.636926,0.638249,0.636978,0.641166,0.637609,0.638667,0.643983,0.642069,0.641121,0.644006,0.644567,0.645366,0.646156,0.647648,0.648130,0.649876,0.648001,0.651155,0.651886,0.651260,0.653817,0.653195,0.655588,0.658066,0.657386,0.658043,0.657431,0.660132,0.660983,0.661772,0.660495,0.664561,0.665503,0.663701,0.664619,0.665684,0.668374,0.668625,0.667248,0.670137,0.671491,0.669408,0.672281,0.676420,0.675892,0.673512,0.675249,0.677192,0.676401,0.678282,0.679802,0.680591,0.682223,0.680239,0.683869,0.684358,0.681643,0.682270,0.682888,0.685516,0.688966,0.688569,0.688257,0.690603,0.690878,0.693088,0.690651,0.693507,0.690976,0.693778,0.695558,0.696775,0.696871,0.697537,0.697071,0.699825,0.701326,0.701488,0.701823,0.702048,0.704433,0.704772,0.702811,0.706289,0.705159,0.707597,0.708045,0.709381,0.709749,0.710342,0.711467,0.711571,0.713882,0.711988,0.714347,0.717014,0.714924,0.717999,0.715959,0.717405,0.720465,0.720147,0.723965,0.722026,0.722034,0.721176,0.723723,0.723843,0.725683,0.729549,0.725753,0.726051,0.728531,0.728295,0.730964,0.733838,0.732339,0.732154,];
Ie=[0.162297,0.159391,0.162831,0.161794,0.164184,0.164706,0.166305,0.165912,0.165029,0.164822,0.165189,0.165125,0.168715,0.165015,0.169773,0.168697,0.169060,0.172074,0.169714,0.169344,0.169910,0.172010,0.173577,0.171981,0.173681,0.171780,0.174392,0.176187,0.173034,0.172210,0.174203,0.173498,0.175625,0.176243,0.175182,0.174879,0.178725,0.175593,0.177061,0.174632,0.176042,0.180057,0.179529,0.177696,0.179875,0.180261,0.181428,0.183351,0.180922,0.178601,0.180453,0.182419,0.184788,0.182757,0.186204,0.187208,0.186650,0.186459,0.185872,0.189257,0.188793,0.187024,0.188688,0.187857,0.189264,0.188977,0.189835,0.191475,0.191092,0.192146,0.191177,0.191983,0.191031,0.191302,0.193094,0.190849,0.191733,0.190528,0.192962,0.199338,0.196476,0.194265,0.192911,0.198851,0.196352,0.200578,0.196255,0.202842,0.199201,0.202425,0.201259,0.201053,0.201620,0.202747,0.201251,0.200675,0.206355,0.205838,0.199444,0.205076,0.203793,0.202268,0.207246,0.210224,0.205404,0.206942,0.205114,0.207374,0.207665,0.207649,0.213519,0.210429,0.210226,0.215097,0.213527,0.211560,0.214066,0.214019,0.216482,0.212610,0.214883,0.213892,0.211533,0.214444,0.213396,0.214236,0.217496,0.216692,0.217683,0.221401,0.218210,0.219348,0.218921,0.221388,0.223638,0.222827,0.225214,0.221808,0.224019,0.223141,0.220978,0.224787,0.226152,0.226480,0.225670,0.226526,0.226186,0.227891,0.229261,0.227309,0.229516,0.227003,0.225570,0.232388,0.233735,0.229088,0.231886,0.233041,0.233884,0.232917,0.236081,0.231952,0.232002,0.233148,0.236094,0.229414,0.234566,0.235379,0.238617,0.239965,0.239242,0.237076,0.238127,0.240483,0.240864,0.239797,0.237987,0.242244,0.244557,0.241332,0.244166,0.243358,0.245840,0.246588,0.242164,0.245247,0.245903,0.241859,0.246696,0.248417,0.244572,0.246920,0.243731,0.249102,0.250442,0.248206,0.251133,0.252693,0.247942,0.254532,0.251959,0.253818,0.252286,0.255125,0.254479,0.254485,0.253605,0.253017,0.257087,0.256205,0.255038,0.255484,0.259662,0.255321,0.260350,0.260636,0.259787,0.262800,0.262532,0.259966,0.262156,0.261079,0.264824,0.263089,0.261912,0.264744,0.266768,0.268112,0.265523,0.265257,0.267791,0.270008,0.269401,0.266930,0.267648,0.268873,0.267510,0.270549,0.268637,0.271158,0.270819,0.271785,0.272388,0.268131,0.276854,0.274643,0.273391,0.273520,0.279220,0.277369,0.277260,0.278663,0.276428,0.279193,0.277638,0.278192,0.280130,0.279313,0.279636,0.280353,0.279079,0.280799,0.281272,0.283595,0.283996,0.284722,0.284123,0.287039,0.283261,0.285185,0.284941,0.287415,0.286373,0.285345,0.287739,0.289468,0.286772,0.287412,0.290832,0.287388,0.289524,0.290123,0.290020,0.293867,0.293268,0.296117,0.290214,0.292118,0.294646,0.295048,0.293002,0.295227,0.294207,0.295337,0.297923,0.297225,0.297104,0.299174,0.295998,0.298459,0.299416,0.299494,0.302639,0.300845,0.302973,0.304615,0.305813,0.303990,0.301780,0.306707,0.304783,0.304643,0.302804,0.306171,0.309721,0.306900,0.306285,0.309149,0.308920,0.311084,0.309283,0.309244,0.312824,0.311060,0.311886,0.315708,0.311681,0.311672,0.312247,0.314540,0.312991,0.315230,0.314171,0.315427,0.314545,0.315196,0.319439,0.317057,0.318032,0.319069,0.317511,0.317789,0.321197,0.320534,0.320748,0.325117,0.320763,0.323582,0.322953,0.324764,0.322132,0.323018,0.325440,0.322191,0.327827,0.324670,0.322796,0.326624,0.326580,0.328618,0.326415,0.327823,0.330086,0.329330,0.327629,0.329824,0.329736,0.330947,0.331341,0.331238,0.330366,0.333488,0.334183,0.333166,0.333120,0.335620,0.333937,0.335866,0.335208,0.336411,0.334202,0.334506,0.338611,0.336634,0.340454,0.339657,0.337363,0.334094,0.339974,0.338159,0.338031,0.341826,0.340800,0.340095,0.339996,0.339996,0.341486,0.343708,0.342640,0.342841,];
p=polyfit(Ia,Ie,3);
xi = 0:0.01:1;
yi = polyval(p,xi);
plot(Ie, Ia,'o','LineWidth',2,'Color',[0 0 0]);
hold on;
plot(yi,xi,'LineWidth',1);
axis([0,1,0,1]);
grid on;
hold off;