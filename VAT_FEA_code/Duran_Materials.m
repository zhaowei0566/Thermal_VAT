switch mat_type
    
    case 'GrEp'
        Mat.E1=155e9;
        Mat.E2=8.07e9;
        Mat.G12=4.55e9;
        Mat.v12=0.22;
        Mat.alpha1= -0.07e-6;
        Mat.alpha2= 30.1e-6;
        
        T01 = [60.7 32.19];
        
        
    case 'EgEp'
        
        
        Mat.E1=41e9;
        Mat.E2=10.04e9;
        Mat.G12=4.3e9;
        Mat.v12=0.28;
        Mat.alpha1= 7e-6;
        Mat.alpha2= 26e-6;
            T01 = [6.710 58.04];  
        
    case 'SgEp'
        Mat.E1=45e9;
        Mat.E2=11e9;
        Mat.G12=4.5e9;
        Mat.v12=0.29;
        Mat.alpha1= 7.1e-6;
        Mat.alpha2= 30.e-6;
        
        T01 = [16.12 54.74];
        
    case 'KeEp'
        Mat.E1=80e9;
        Mat.E2=5.5e9;
        Mat.G12=2.2e9;
        Mat.v12=0.34;
        Mat.alpha1= -2e-6;
        Mat.alpha2= 60e-6;
        
       T01 = [66.05 11.73];
        
        
    case 'CaEp'
        Mat.E1=147e9;
        Mat.E2=10.3e9;
        Mat.G12=7e9;
        Mat.v12=0.27;
        Mat.alpha1= -0.9e-6;
        Mat.alpha2= 27e-6;
        
          T01 = [69 -5.705];
        
    case 'CaPe'
        Mat.E1=138e9;
        Mat.E2=8.7e9;
        Mat.G12=5e9;
        Mat.v12=0.28;
        Mat.alpha1= -0.2E6;
        Mat.alpha2= 24e-6;
        
           T01 = [63.07 29.5];
        
    case 'CaPo'
        Mat.E1=216e9;
        Mat.E2=5e9;
        Mat.G12=4.5e9;
        Mat.v12=0.25;
        Mat.alpha1= -0.0e-6;
        Mat.alpha2= 25e-6;
        
        T01 = [56.3 36.68];
        
    case 'BoEp'
        
        Mat.E1=201e9;
        Mat.E2=21.7e9;
        Mat.G12=5.4e9;
        Mat.v12=0.17;
        Mat.alpha1= 6.1e-6;
        Mat.alpha2= 30.e-6;
        
        T01 = [-6.57 63.28];
        
        
end