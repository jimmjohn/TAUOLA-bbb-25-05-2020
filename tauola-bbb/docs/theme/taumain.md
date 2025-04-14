#  taumain



In this section, we'll explain how to call the main Tauola subroutines and use them to simulate tau decays. To use the Tauola library, certain COMMON blocks must be initialized to specify necessary parameters.

 ***We need to move these constants and decay initialization to a separate file, as the code that the user modifies should typically contain only the information relevant to their specific application.***

The following initialization's has to be done before calling the `DEXAY`/`DEKAY` function from the `Tauola` library.

The `IDFF` in `IDFC` common block should be initiated to 15, if we want to run the decay for \f$\tau^+\f$ decay and and -15 for the \f$\tau^-\f$ decay.

#### Mass & Particle ID convention in Tauola(DCDMAS/INIMAS)

|      Particle      |    Mass     |   ID   | Variable |
| :----------------: | :---------: | :----: | :------: |
|    \f$m_\tau\f$    |    1.777    |        |  AMTAU   |
| \f$m_{\nu_\tau}\f$ |     0.0     |        |  AMNUTA  |
|     \f$m_e\f$      |  0.000511   |   11   |   AMEL   |
|  \f$m_{\nu_e}\f$   |     0.0     |   12   |  AMNUE   |
|    \f$m_\mu\f$     |  0.105659   |   13   |   AMMU   |
|  \f$m_{\pi^0}\f$   |  0.134976   | 2, 111 |  AMPIZ   |
|   \f$m_{\pi}\f$    |  0.139570   |   1    |   AMPI   |
|    \f$m_\rho\f$    |   0.77590   |        |   AMRO   |
|   \f$m_{a_1}\f$    |    1.251    |        |   AMA1   |
|     \f$m_K\f$      |  0.493677   |   3    |   AMK    |
|   \f$m_{K^0}\f$    |  0.497672   |   4    |   AMKZ   |
|   \f$m_{K^*}\f$    |   0.89166   |  313   |  AMKST   |
|    \f$m_{GM}\f$    |   0.0001    |   8    |          |
|    \f$m_\eta\f$    |   0.5488    | 9, 221 |          |
| \f$m_{\nu_\mu}\f$  |     0.0     |   14   |  AMNUMU  |
|  \f$m_{\gamma}\f$  |     0.0     |   22   |          |
|   \f$m_\omega\f$   |   0.7826    |  223   |          |
|    \f$m_\phi\f$    |    1.019    |  333   |          |
|  \f$m_{\rho^0}\f$  |   0.7755    |  113   |          |
|  \f$m_{\eta_p}\f$  |    0.958    |  331   |          |
|     \f$m_p\f$      | 0.938272046 |  2212  |          |
|   \f$m_{A^0}\f$    |    1.45     | 10211  |          |
|   \f$m_{B_1}\f$    |    1.235    | 10213  |          |
|   \f$m_{K_s}\f$    |  0.497614   |  310   |          |
|  \f$m_\lambda\f$   |   1.15683   |  3122  |          |
|   \f$m_{a_0}\f$    |    0.98     | 10111  |          |
|   \f$m_{f_0}\f$    |    0.98     | 10221  |          |



***Doubt -*** 

1. ​	The PDG ID of \f$K^0\f$ is 311 and in the code its assigned a mass of \f$\pi^0\f$(line number 1388 in taumain). 



#### TAUDCDsize.in

Following values in the TAUDCDsize.in file is used to define the number of different channels in different modes.

| Variable | Value | Description                                        |
| :------: | :---: | -------------------------------------------------- |
| `NMODE`  |  196  | Maximum possible number of decay channels          |
|  `NLT`   |   2   | Number of leptonic decay channels                  |
|  `NM1`   |  40   | Maximum number of 1 scalar or *anomalous* channels |
|  `NM2`   |  71   | Maximum number of 2 scalar or *anomalous* channels |
|  `NM3`   |  19   | Maximum number of 3 scalar channels                |
|  `NM4`   |  32   | Maximum number of 4 scalar channels                |
|  `NM5`   |  21   | Maximum number of 5 scalar channels                |
|  `NM6`   |  13   | Maximum number of multiple scalar channels         |

The definition of different channels follows the following order

1. Leptonic

2. 4 - scalar channels

3. 5 - scalar channels

4. multiple scalar channels

5. 3 - scalar channels

6. 2 - scalar channels

7. 1 - scalar or *anomalous* channels

   



#### INITDK

This routine initializes the decay channels and the branching ratios.

###### The channel number is saved to `JLIST` and the corresponding branching ratio is saved to `GAMPRT` array. Also the IDFFIN(J, JNPI) is saved with the particle ID's for the channel (nchannel-2) from J=0 to J=JMAX. The maximum particles in a channel is set to 9.

##### The leptonic channel Branching ratio's

1. 
   \f$\tau \to e^-\f$ , GAMPRT = 0.1800(CLEO default)  , GAMPRT = 0.178651(BaBar)
2. \f$\tau \to \mu^-\f$ , GAMPRT = 0.1751(CLEO default)  , GAMPRT = 0.173551(BaBar)
   

##### Four Scalar Channels

1. \f$\tau^- \to  2\pi^- + \pi^+ + \pi^0\f$ , GAMPRT = 0.0450, GAMPRT = 0.043654(BaBar)
2. \f$\tau^- \to  3\pi^0 + \pi^-\f$ , GAMPRT = 0.0100,  GAMPRT = 0.012619(BaBar)
3. \f$\tau^- \to  \nu_e + 2e^- + e^+\f$ , GAMPRT = 0.0100 * 0
4. \f$\tau^- \to  \nu_\mu + 2\mu^- + \mu^+\f$ , GAMPRT = 0.0100 * 0
5. \f$\tau^- \to  \nu_e + e^- + \mu^- + \mu^+\f$ , GAMPRT = 0.0100 * 0
6. \f$\tau^- \to  \nu_\mu + \mu^- + e^- + e^+\f$ , GAMPRT = 0.0100 * 0
7. \f$\tau^- \to  K^- + 3\pi^0\f$ , GAMPRT = 0.0100 * 0
8. \f$\tau^- \to  2\pi^0 + \eta + K^-\f$ , GAMPRT = 0.0100 * 0
9. \f$\tau^- \to  2\pi^0 + K^0 + \pi^-\f$ , GAMPRT = 0.0100 * 0
10. \f$\tau^- \to  \pi^0 + K^0 + \eta + \pi^-\f$ , GAMPRT = 0.0100 * 0
11. \f$\tau^- \to  \pi^0 + \pi^- + \pi^+ + K^-\f$ , GAMPRT = 0.0100 * 0
12. \f$\tau^- \to  K^0 + \pi^- + \pi^+ + \pi^-\f$ , GAMPRT = 0.0100 * 0
13. \f$\tau^- \to  2\pi^0 + \eta + \pi^-\f$ , GAMPRT = 0.0100 * 0
14. \f$\tau^- \to  K^0 + \bar{K^0} + \eta + \pi^-\f$ , GAMPRT = 0.0100 * 0
15. \f$\tau^- \to  K^0 + \bar{K^0} + \pi^0 + \pi^-\f$ , GAMPRT = 0.0100 * 0
16. \f$\tau^- \to  K^0 + \bar{K^0} + K^0 + \pi^-\f$ , GAMPRT = 0.0100 * 0
17. \f$\tau^- \to  K^0 + 2\pi^0 + K^-\f$ , GAMPRT = 0.0100 * 0
18. \f$\tau^- \to  K^0 + \bar{K^0} + \pi^0 + K^-\f$ , GAMPRT = 0.0100 * 0
19. \f$\tau^- \to  \pi^0 + K^0 + \eta + K^-\f$ , GAMPRT = 0.0100 * 0
20. \f$\tau^- \to  \pi^- + \pi^+ + \pi^- + \eta\f$ , GAMPRT = 0.0100 * 0
21. \f$\tau^- \to  \pi^- + K^+ + K^- + \pi^0\f$ , GAMPRT = 0.0100 * 0
22. \f$\tau^- \to  K^- + K^+ + K^- + \pi^0\f$ , GAMPRT = 0.0100 * 0
23. \f$\tau^- \to  K^- + K^+ + K^- + K^0\f$ , GAMPRT = 0.0100 * 0
24. \f$\tau^- \to  K^- + \pi^+ + \pi^- + K^0\f$ , GAMPRT = 0.0100 * 0
25. \f$\tau^- \to  K^- + K^+ + \pi^- + K^0\f$ , GAMPRT = 0.0100 * 0
26. \f$\tau^- \to  \pi^- + \pi^+ + \pi^- + \omega\f$ , GAMPRT = 0.0100 * 0
        

There are another 6 more placeholder for including 6 more decay channels so total will be 32 (NM4=32)\\


The ID of the particle decayed will be saved in `IDFFIN`(particleNumber, Channel)

\f$ \tau^- \to  2\pi^- + \pi^+ + \pi^0\f$

- IDFFIN(1,1) = -1

- IDFFIN(2,1) = -1

- IDFFIN(3,1) = 1

- IDFFIN(4,1) = 2

- IDFFIN(5,1) = 0

- IDFFIN(6,1) = 0

- IDFFIN(7,1) = 0
- IDFFIN(8,1) = 0
- IDFFIN(9,1) = 0

\f$\tau^- \to  3\pi^0 + \pi^-\f$

- IDFFIN(1,1) = 2
- IDFFIN(2,1) = 2
- IDFFIN(3,1) = 2
- IDFFIN(4,1) = -1
- IDFFIN(5,1) = 0
- IDFFIN(6,1) = 0
- IDFFIN(7,1) = 0
- IDFFIN(8,1) = 0
- IDFFIN(9,1) = 0

\f$\tau^- \to  \nu_e + 2e^- + e^+\f$

- IDFFIN(1,1) = 12
- IDFFIN(2,1) = -11
- IDFFIN(3,1) = -11
- IDFFIN(4,1) = 11
- IDFFIN(5,1) = 0
- IDFFIN(6,1) = 0
- IDFFIN(7,1) = 0
- IDFFIN(8,1) = 0
- IDFFIN(9,1) = 0



##### Five scalar channels

1. \f$\tau^- \to 2\pi^- + \pi^+ + 2\pi^0 \;old\f$ , GAMPRT=0.0009, GAMPRT = 0.005011(BaBar)
2. \f$\tau^- \to a_1 \to \rho + \omega\f$, GAMPRT=0.00
3. \f$\tau^- \to benchmark \;current\f$, GAMPRT=0.0
4. \f$\tau^- \to 2\pi^- + \pi^+ + 2\pi^0  \;app08\f$  , GAMPRT=0.00
5. \f$\tau^- \to \pi^- + 4\pi^0 \;app08\f$     ,GAMPRT=0.00
6. \f$\tau^- \to 3\pi^- + 2\pi^+ \;app08\f$     ,GAMPRT=0.00
7. \f$\tau^- \to 2\pi^- + 2\pi^+ + K^-\f$       ,GAMPRT=0.001*0
8. \f$\tau^- \to 2\pi^- + \pi^+ + \pi^0 + K^0\f$
9. \f$\tau^- \to \pi^- +  4\pi^0 \;old\f$



##### Multiple scalar channels

1. \f$\tau^- \to 3 \pi^- + 2\pi^+\f$,  GAMPRT=0.0004, GAMPRT = 0.000789(BaBar)
2. \f$\tau^- \to 3 \pi^- + 2\pi^+ + \pi^0\f$,  GAMPRT=0.0003, GAMPRT = 0.000183(BaBar)
3. \f$\tau^- \to 2 \pi^- + \pi^+ + 3\pi^0\f$,  GAMPRT=0.0005, GAMPRT = 0.000251(BaBar)
4. \f$\tau^- \to 3 \pi^- + 2\pi^+ + 2\pi^0\f$,  GAMPRT=0.0005*0
5. \f$\tau^- \to 4 \pi^- + 3\pi^+\f$,  GAMPRT=0.0005*0
6. \f$\tau^- \to 4 \pi^- + 3\pi^+ + \pi^0\f$,  GAMPRT=0.0005*0
7. \f$\tau^- \to 2 \pi^- + 2\pi^+ + K^- + \pi^0\f$,  GAMPRT=0.0005*0



##### 3 - scalar channels

1. \f$\tau^- \to K^- + \pi^- + K^+\f$,  GAMPRT=0.0015,  GAMPRT =0.00159 (BaBar)
2. \f$\tau^- \to K^0 + \pi^- + K^0_B\f$,  GAMPRT=0.0015,  GAMPRT = 0.001672(BaBar)
3. \f$\tau^- \to K^- + \pi^0 + K^0\f$,  GAMPRT=0.0015,  GAMPRT = 0.001536(BaBar)
4. \f$\tau^- \to \pi^0 + \pi^0 + K^-\f$,  GAMPRT=0.0005,  GAMPRT = 0.00068(BaBar)
5. \f$\tau^- \to K^- + \pi^- + \pi^+\f$,  GAMPRT=0.0050,  GAMPRT = 0.003009(BaBar)
6. \f$\tau^- \to \pi^- + K^0_B + \pi^0\f$,  GAMPRT=0.0055,  GAMPRT = 0.003767(BaBar)
7. \f$\tau^- \to \eta + \pi^- + \pi^0\f$,  GAMPRT=0.0017,  GAMPRT = 0.00183(BaBar)
8. \f$\tau^- \to \pi^0 + \pi^0 + \gamma \f$,  GAMPRT=0.0013,  GAMPRT = 0.000802(BaBar)
9. \f$\tau^- \to \pi^0 + \pi^0 + \pi^-\f$,  GAMPRT=0.1790/2,  GAMPRT = 0.091783(BaBar)
10. \f$\tau^- \to \pi^- + \pi^- + \pi^+\f$,  GAMPRT=0.1790/2,  GAMPRT = 0.091783(BaBar)
11. \f$\tau^- \to K^- + K^- + K^+\f$,  GAMPRT=0.0010 *0
12. \f$\tau^- \to K^- + K^0 + K^0\f$,  GAMPRT=0.0010*0
13. \f$\tau^- \to K^- + \eta + \pi^0\f$,  GAMPRT=0.0010*0
14. \f$\tau^- \to K^0 + \eta + \pi^-\f$,  GAMPRT=0.0010*0
15. \f$\tau^- \to K^- + K^0 + \rho^0\f$,  GAMPRT=0.0010*0
16. \f$\tau^- \to \pi^- + \phi + \pi^0\f$,  GAMPRT=0.0010*0
17. \f$\tau^- \to K^- + \phi + \pi^0\f$,  GAMPRT=0.0010*0
18. \f$\tau^- \to K^0 + \eta + K^-\f$,  GAMPRT=0.0010*0



##### 2 - scalar channels

1. \f$\tau^- \to \pi^- + \pi^0\f$,  GAMPRT=0.2515,  GAMPRT = 0.253754(BaBar)
2. \f$\tau^- \to \pi^- + K^0\f$,  GAMPRT=0.0134\*0.6666,  GAMPRT = 0.013641\*0.6666(BaBar)
3. \f$\tau^- \to K^- + \pi^0\f$,  GAMPRT=0.0134\*0.3334,  GAMPRT = 0.013641\*0.3334(BaBar)
4. \f$\tau^- \to K^- + K^0\f$,  GAMPRT=0.0010,  GAMPRT = 0.001651(BaBar)

##### 1 - scalar channels

1. \f$\tau^- \to \pi^- \f$, GAMPRT=0.1110,  GAMPRT = 0.110841(BaBar)
2. \f$\tau^- \to K^- \f$, GAMPRT=0.0071,  GAMPRT = 0.006946(BaBar)



The \f$K^0\f$ swapped between \f$K_{L}\f$ and \f$K_{S}\f$ 50 % of the time, with \f$K_L\f$ particle ID 130 and \f$K_S\f$ particle ID 310.



#### INIPHY

It defines certain constants like ALFINV=1/\f$\alpha\f$, ALFPI = \f$\pi \alpha\f$ and it sets the minimum threshold of the energy scale(XK0), in this case we set it at `XK0`= 0.001



#### INISAMPL

Initialize various parameters used in optimization of phase space generation. The probabilities for different resonances is passed through the `SAMPL2`, `SAMPL3`, `SAMPL4`, and `SAMPL5` common block for different scalar modes.

- PROB1 - ?
- PROB2 - ?
- PROB3 - ?

##### INSAMPL2 is the initialization for 2 scalar, there are 71 modes with 2 scalar

AM2(1), AM2(2)....., AM(71) = \f$m_{K^*}\f$

AM3(1), AM3(2)....., AM(71) = \f$m_\rho\f$

GAM2(1), GAM(2),...., GAM(71) = \f$\Gamma_{K^*}\f$

GAM3(1), GAM3(2),...., GAM3(71) = \f$\Gamma_{\rho}\f$

- For first channel PROB1=0, PROB2=0.
- For second and third channels PROB1=0, PROB2=1.
- For other channels PROB1=1, PROB2=0.

Note: This PROB1, PROB2 is in SAMPL2

##### INSAMPL3 is the  initialization for 3 scalar, there are 19 modes with 3 scalar

\f$m_{\rho_p}\f$ = 1.1

\f$\Gamma_{\rho_p}\f$ = 0.36

\f$m_\omega\f$ = 0.782

\f$\Gamma_\omega\f$ = 0.0084

Probabilities of different channels

1. PROB1=0.5, PROB2=0.5, AMRX=1.57, GMRX=0.9, AMRA=\f$m_{K^*}\f$, AMRB=\f$m_\rho\f$
2. PROB1=0.5, PROB2=0.5, AMRX=1.57, GMRX=0.9, AMRA=\f$m_{K^*}\f$, AMRB=\f$m_\rho\f$
3. PROB1=0.5, PROB2=0.5, AMRX=1.27, GMRX=0.3, AMRA=\f$m_{K^*}\f$, AMRB=\f$m_{K^*}\f$
4. PROB1=0.5, PROB2=0.5, AMRX=1.27, GMRX=0.3, AMRA=\f$m_{K^*}\f$, AMRB=\f$m_{K^*}\f$
5. PROB1=0.5, PROB2=0.5, AMRX=1.27, GMRX=0.3, AMRA=\f$m_{K^*}\f$, AMRB=\f$m_\rho\f$
6. PROB1=0.4, PROB2=0.4, AMRX=1.27, GMRX=0.3, AMRA=\f$m_\rho\f$, AMRB=\f$m_{K^*}\f$
7. PROB1=0., PROB2=1., AMRX=1.27, GMRX=0.9, AMRA=\f$m_\rho\f$, AMRB=\f$m_\rho\f$
8. PROB1=0., PROB2=1., AMRX=\f$m_{\rho_p}\f$, AMRA=\f$m_\omega\f$, AMRB=\f$m_\rho\f$
9. PROB1=0.5, PROB2=0.5, AMRX=\f$m_{a1}\f$, AMRA=\f$m_\rho\f$, AMRB=\f$m_\rho\f$   
10. PROB1=0.5, PROB2=0.5, AMRX=\f$m_{a1}\f$, AMRA=\f$m_\rho\f$, AMRB=\f$m_\rho\f$

For all others

​    PROB1=0.0, PROB2=0.0, AMRX=\f$m_{a1}\f$, AMRA=\f$m_\rho\f$, AMRB=\f$m_\rho\f$


##### INSAMPL4 is the  initialization for 4 scalar, there are 32 modes with 4 scalar

Probabilities of different modes

1- PROB1=0.35, PROB2=0.35, AMRX=1.2, GAMRX=0.46, AMRA=\f$m_\omega\f$

2 - PROB1=0.0, PROB2=0.0, AMRX=1.4, GAMRX=0.6, AMRA=\f$m_\omega\f$

3 -- 12. PROB1=0.0, PROB2=0.0, AMRX=1.4, GAMRX=0.6, AMRA=\f$m_\omega\f$

For all others

​	PROB1=0.0, PROB2=0.0, AMRX=\f$m_{a1}\f$, AMRA=\f$m_\rho\f$


##### INSAMPL5 is the  initialization for 5 scalar, there are 21 modes with 5 scalar

Probabilities of all modes

PROBa2 = 0.7

PROBOM = 0.7

ama2 = 1.260

gama2 = 0.4

\f$m_\omega\f$ = 0.78257

\f$\Gamma_\omega\f$ = 0.7

For 1\f$^{st}\f$ and 2\f$^{nd}\f$ modes, \f$\Gamma_\omega\f$ = 0.00844

#### TAUFIL

This subroutine is responsible for simulating tau lepton production and storing its momentum and identity in a LUND common block. 

XPB1, XPB2 - Beam momenta. 
AQF1, AQF2 - \f$\tau^+\f$ and \f$\tau^-\f$ four momenta.



###### Some of the variables like nhep is not initialized. This will lead to uncertain results when we call FILHEP 


|                        Read Next |
|---------------------------------:|
| [Extensions](docs/extensions.md) |

</div>
