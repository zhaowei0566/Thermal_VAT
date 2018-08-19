function fiber_ply = VAT_fiber(T01,number,flag)


fiber_ply = zeros(2*number,2);

t0  = T01(1);
t1  = T01(2);

half_laminates = zeros(number,2);


for ii = 1:number
    
    
    half_laminates(ii,:) = [t0 t1]*(-1)^(ii-1);
    
end


switch flag
    
    
    case 'ASYM'
        
        
        fiber_ply (1:number,:) =   half_laminates;
        
        fiber_ply (1+number:end,:) =   half_laminates;
        
    case 'SYM'
        
        
                fiber_ply (1:number,:) =   half_laminates;
        
        fiber_ply (1+number:end,:) =   flipud(half_laminates);
end

