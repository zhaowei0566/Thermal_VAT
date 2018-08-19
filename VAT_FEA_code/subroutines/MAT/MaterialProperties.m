function Mat=MaterialProperties(MaterialName)

switch MaterialName
    
    case 'GraphitePolymer'
        
        Mat.E1=155*10^9;
        Mat.E2=12.1*10^9;
        Mat.E3=12.1*10^9;
        
        Mat.v12=0.248;
        Mat.v13=0.248;
        Mat.v23=0.458;
        
        Mat.G12=4.40*10^9;
        Mat.G13=4.40*10^9;
        Mat.G23=3.20*10^9;
        
        Mat.alpha1=-0.018*10^(-6);
        Mat.alpha2=24.3*10^(-6);
        Mat.alpha3=24.3*10^(-6);
        
        Mat.beta1=146*10^(-6);
        Mat.beta2=4770*10^(-6);
        Mat.beta3=4770*10^(-6);
        
        % In-plane yielding stress in 1- and 2- axis
        Mat.stress1C=-1250*10^6;
        Mat.stress1T=1500*10^6;
        Mat.stress2C=-200*10^6;
        Mat.stress2T=50*10^6;
        
        Mat.tau12F=100*10^6;
        
        Mat.layerthickness=150*10^(-6);
        
        Mat.density=1800;
        
        
    case 'GlassPolymer'
        
        Mat.E1=50*10^9;
        Mat.E2=15.2*10^9;
        Mat.E3=15.2*10^9;
        
        Mat.nu12=0.254;
        Mat.nu13=0.254;
        Mat.nu23=0.428;
        
        Mat.G12=4.70*10^9;
        Mat.G13=4.70*10^9;
        Mat.G23=3.28*10^9;
        
        Mat.alpha1=6.34*10^(-6);
        Mat.alpha2=23.3*10^(-6);
        Mat.alpha3=23.3*10^(-6);
        
        Mat.beta1=434*10^(-6);
        Mat.beta2=6320*10^(-6);
        Mat.beta3=6320*10^(-6);
        
        
        % In-plane yielding stress in 1- and 2- axis
        Mat.stress1C=-600*10^6;
        Mat.stress1T=1000*10^6;
        Mat.stress2C=-120*10^6;
        Mat.stress2T=30*10^6;
        Mat.tau12F=70*10^6;
        
        
    case 'Aluminum'
        
        Mat.E1=72.4*10^9;
        Mat.E2=72.4*10^9;
        Mat.E3=72.4*10^9;
        
        Mat.nu12=0.3;
        Mat.nu13=0.3;
        Mat.nu23=0.3;
        
        Mat.G12=Mat.E1/2/(1+Mat.nu12);
        Mat.G13=Mat.E1/2/(1+Mat.nu12);
        Mat.G23=Mat.E1/2/(1+Mat.nu12);
        
        Mat.alpha1=22.5*10^(-6);
        Mat.alpha2=22.5*10^(-6);
        Mat.alpha3=22.5*10^(-6);
        
        Mat.beta1=0;
        Mat.beta2=0;
        Mat.beta3=0;
        
    case 'ThinWalledBeam'
        % Chandra, Stemple and Chopra Thin-walled composite beam
        Mat.E1=20.59e6;
        Mat.E2=1.42e6;
        
        Mat.G12=0.89e6;
        
        Mat.nu12=0.42;
        
       Mat.alpha1=22.5*10^(-6);
        Mat.alpha2=22.5*10^(-6);
        Mat.alpha3=22.5*10^(-6);
        
        Mat.beta1=0;
        Mat.beta2=0;
        Mat.beta3=0;
        
end

% Tsai-Wu Constants
% Mat.F1=(1/Mat.stress1T+1/Mat.stress1C);
% Mat.F2=(1/Mat.stress2T+1/Mat.stress2C);
% Mat.F66=(1/Mat.tau12F)^2;
% Mat.F11=-1/(Mat.stress1T*Mat.stress1C);
% Mat.F22=-1/(Mat.stress2T*Mat.stress2C);
