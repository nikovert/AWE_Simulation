classdef AWE_3DOF < DynSys
    properties
        ENVMT           % - environmental variable
        AIRCRAFT        % - aircarft variables
        T               % - tether properties 
        base_windspeed  % - wind speed at 6m
        Ft_set          % - set tether force (used for simpliefied tether dynamics)
        v_ro_set        % - set reelout speed
        
        LongLatState = false   % - if true, first two states are expected to be long, lat
        
        F_T_max         % maximum tether force before rupture
        max_tether_diff_dot = 0.005;  % - max reel-out/in speed
        d_wind_max = 4; % max wind turbulence velocity

        % To be replaced later on
        v_w_O
        
        % The maximum angle
        alpha_max = pi/18
        mu_max = 1
        alpha_options = 5; %WARNING MUST UPDATE OPTCTRL IF CHANGED HERE
        mu_options = 9;
        
        % remove tether forces for testing purposes
        skipTether = false
        
        % To speed up optCtrl
        optCtrlMode = 'normal'
        alpha
        mu
        alphaMax
        muMax
        alphaMin
        muMin
        
        % direction for following the curve
        curve_direction
        Lem
        
        % Normalization
        h0  % distance 
        v0  % speed
        a0  % angle
        
        % To speed up dynamics
        F_rest
        F_tether_drag_Abar
        pos_W
        long_dot
        lat_dot
        h_tau_dot
        s_dot
        sigma_dot
        final_seg_diff_dot
        
        ignoreTetherDiff % if true, we assume that teh tetherdiff stayes constant
    end % end properties

    methods
        function obj = AWE_3DOF(x, base_windspeed, ENVMT, AIRCRAFT, T)
        % obj = AWE3DF(x)
        %
        %   x is expected to already be normalised
        %
        % ENVMT     - environmental variable
        % AIRCRAFT  - aircarft variables
        % T         - tether properties 
        %
        % STATES
        %     long    - longitude defined in the wind reference frame 
        %     lat     - latitude defined in the wind reference frame
        %     h_tau   - radial distance to the origin in the wind reference frame
        %     va      - apparent wind speed in the A_bar (rotated aerodynamic frame) frame
        %     chi_a   - course angle in the A_bar frame
        %     gamma_a - path angle in the A_bar frame
        %     l_s     - tether segment length
        %
        % INPUTS
        %     alpha_a - aerodynamic bank angle defined in the A_bar frame
        %     mu_a    - aerodynamic bank angle
        %     
        %     The sideslip angle (beta_a) is set to zero for the 3DOF model
        %
        %

        if nargin < 5
            %% Tether parameter
            T.c0 =  614600 ;
            T.d0 = 473;
            T.rho_t = 0.518/170;
            T.d_tether = 0.004;
            T.CD_tether = 1.2;
            T.np = 5;
            T.CD_tether = 1.2;
            T.E = 5.3e9;
            T.eps = 0.02;
            T.springconstant_un = T.E*pi*T.d_tether^2/(4*0.02);
            T.d_tether = 2e-3;
            T.l0 = (x(3)+1)/(T.np+1);
        end

        if nargin < 4
            %% Aircraft parameters
            AIRCRAFT.AP_Ts= 0.0100;
            AIRCRAFT.delay_g2k= 0.1000;
            AIRCRAFT.delay_k2g= 0.1000;
            AIRCRAFT.S_ref= 3;
            AIRCRAFT.S_wing= 3;
            AIRCRAFT.b= 5.5000;
            AIRCRAFT.c= 0.5500;
            AIRCRAFT.mass= 36.8000;
            AIRCRAFT.J= [25  , 0 , -0.47; ...
                           0 , 32,  0   ; ...
                        -0.47, 0 , 56  ];
            AIRCRAFT.Jx= 25;
            AIRCRAFT.Jy= 32;
            AIRCRAFT.Jz= 56;
            AIRCRAFT.Jxz= 0.4700;
            AIRCRAFT.Jinv = inv( AIRCRAFT.J ); 
            AIRCRAFT.d_tether= 0.0025;
            AIRCRAFT.rho_tether= 0.0046;
            AIRCRAFT.Cd_tether= 1.2000;
            AIRCRAFT.rho_air= 1.2250;
            AIRCRAFT.Cx_0_0= -0.0293;
            AIRCRAFT.Cx_0_alpha= 0.4784;
            AIRCRAFT.Cx_0_alpha2= 2.5549;
            AIRCRAFT.Cx_q_0= -0.6029;
            AIRCRAFT.Cx_q_alpha= 4.4124;
            AIRCRAFT.Cx_q_alpha2= 0;
            AIRCRAFT.Cx_deltaE_0= -0.0106;
            AIRCRAFT.Cx_deltaE_alpha= 0.1115;
            AIRCRAFT.Cx_deltaE_alpha2= 0;
            AIRCRAFT.Cy_beta_0= -0.1855;
            AIRCRAFT.Cy_beta_alpha= -0.0299;
            AIRCRAFT.Cy_beta_alpha2= 0.0936;
            AIRCRAFT.Cy_p_0= -0.1022;
            AIRCRAFT.Cy_p_alpha= -0.0140;
            AIRCRAFT.Cy_p_alpha2= 0.0496;
            AIRCRAFT.Cy_r_0= 0.1694;
            AIRCRAFT.Cy_r_alpha= 0.1368;
            AIRCRAFT.Cy_r_alpha2= 0;
            AIRCRAFT.Cy_deltaA_0= -0.0514;
            AIRCRAFT.Cy_deltaA_alpha= -0.0024;
            AIRCRAFT.Cy_deltaA_alpha2= 0.0579;
            AIRCRAFT.Cy_deltaR_0= 0.1032;
            AIRCRAFT.Cy_deltaR_alpha= 0.0268;
            AIRCRAFT.Cy_deltaR_alpha2= -0.1036;
            AIRCRAFT.Cz_0_0= -0.5526;
            AIRCRAFT.Cz_0_alpha= -5.0676;
            AIRCRAFT.Cz_0_alpha2= 5.7736;
            AIRCRAFT.Cz_q_0= -7.5560;
            AIRCRAFT.Cz_q_alpha= 0.1251;
            AIRCRAFT.Cz_q_alpha2= 6.1486;
            AIRCRAFT.Cz_deltaE_0= -0.3150;
            AIRCRAFT.Cz_deltaE_alpha= -0.0013;
            AIRCRAFT.Cz_deltaE_alpha2= 0.2923;
            AIRCRAFT.Cm_0_0= -0.0307;
            AIRCRAFT.Cm_0_alpha= -0.6027;
            AIRCRAFT.Cm_0_alpha2= 0;
            AIRCRAFT.Cm_q_0= -11.3022;
            AIRCRAFT.Cm_q_alpha= -0.0026;
            AIRCRAFT.Cm_q_alpha2= 5.2885;
            AIRCRAFT.Cm_deltaE_0= -1.0427;
            AIRCRAFT.Cm_deltaE_alpha= -0.0061;
            AIRCRAFT.Cm_deltaE_alpha2= 0.9974;
            AIRCRAFT.Cl_beta_0= -0.0630;
            AIRCRAFT.Cl_beta_alpha= -3.0000e-04;
            AIRCRAFT.Cl_beta_alpha2= 0.0312;
            AIRCRAFT.Cl_p_0= -0.5632;
            AIRCRAFT.Cl_p_alpha= -0.0247;
            AIRCRAFT.Cl_p_alpha2= 0.2813;
            AIRCRAFT.Cl_r_0= 0.1811;
            AIRCRAFT.Cl_r_alpha= 0.6448;
            AIRCRAFT.Cl_r_alpha2= 0;
            AIRCRAFT.Cl_deltaA_0= -0.2489;
            AIRCRAFT.Cl_deltaA_alpha= -0.0087;
            AIRCRAFT.Cl_deltaA_alpha2= 0.2383;
            AIRCRAFT.Cl_deltaR_0= 0.0044;
            AIRCRAFT.Cl_deltaR_alpha= -0.0013;
            AIRCRAFT.Cl_deltaR_alpha2= 0;
            AIRCRAFT.Cn_beta_0= 0.0577;
            AIRCRAFT.Cn_beta_alpha= -0.0849;
            AIRCRAFT.Cn_beta_alpha2= 0;
            AIRCRAFT.Cn_p_0= -0.0565;
            AIRCRAFT.Cn_p_alpha= -0.9137;
            AIRCRAFT.Cn_p_alpha2= 0;
            AIRCRAFT.Cn_r_0= -0.0553;
            AIRCRAFT.Cn_r_alpha= 0.0290;
            AIRCRAFT.Cn_r_alpha2= 0.0257;
            AIRCRAFT.Cn_deltaA_0= 0.0190;
            AIRCRAFT.Cn_deltaA_alpha= -0.1147;
            AIRCRAFT.Cn_deltaA_alpha2= 0;
            AIRCRAFT.Cn_deltaR_0= -0.0404;
            AIRCRAFT.Cn_deltaR_alpha= -0.0117;
            AIRCRAFT.Cn_deltaR_alpha2= 0.0409;
            AIRCRAFT.CL0= 0.4687;
            AIRCRAFT.CL_alpha= 4.5619;
            AIRCRAFT.g= 9.8100;
            AIRCRAFT.sampling_rate_fcs= 100;
            AIRCRAFT.k_motor= 80;
            AIRCRAFT.S_prop= 0.2027;
            AIRCRAFT.C_prop= 1;
            AIRCRAFT.maxBandwidth= 35;
        end

        if nargin < 3
            %% Environment struct
            ENVMT.g= 9.8066;
            ENVMT.rhos= 1.2250;
            ENVMT.windDirection_rad= 3.1416;
        end
        
        if nargin < 2
            base_windspeed = 9;
        end
        obj.T = T;
        obj.AIRCRAFT = AIRCRAFT;
        obj.ENVMT = ENVMT;
        obj.base_windspeed = base_windspeed;
        
        %% Process initial state
        obj.x = x;
        obj.xhist = x;
        
        obj.nx = length(x);
        obj.nu = 2;
        end % end constructor
        
        function reset_dynamics(obj)
            obj.F_rest = [];
            obj.F_tether_drag_Abar = [];
            obj.pos_W = [];
            obj.long_dot = [];
            obj.lat_dot = [];
            obj.h_tau_dot = [];
            obj.ls_diff_dot = [];
        end
    end % end methods
end % end class