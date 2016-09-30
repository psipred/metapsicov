#define IPERGRP (27)

#define WINR (7)
#define WINL (-WINR)

#define CWINR (5)
#define CWINL (-CWINR)

#define NUM_IN	(2*(WINR-WINL+1)*IPERGRP+(CWINR-CWINL+1)*IPERGRP+(CWINR-CWINL+1)*(CWINR-CWINL+1)*6+16+21+8)	/* number of input units */
#define NUM_HID (160)			/* number of hidden units */
#define NUM_HID2 (160)			/* number of hidden units */
#define NUM_OUT (2) 			/* number of output units */

#define TOTAL		(NUM_IN + NUM_HID + NUM_HID2 + NUM_OUT)

#define IALPHA		(0.9)		/* Initial smoothing factor */
#define ILRATE		(0.001)         /* Initial learning rate */
#define RWMAX		(0.03)		/* Max random weight value */
#define SAMPINTV	(1)		/* Sampling interval */
#define SAMPSIZE	(0)		/* Samples per weight update */

#define CPU_MAX		(10*3600.0)	/* Maximum CPU time (secs) for learning */

#define EPOCHMODE
