function H_s = val_Hellinger(Resp_r, Resp_s, npoints)
% computes Hellinger distance between two noramlly distributed samples

% INPUT:
% Resp_r = reference response
% Resp_s = sampled response
% npoints = number of integration points 

% OUTPUT:
% H_s = Hellinger distance
    

    %initialize variables
    H_s = zeros(size(Resp_r,1),1);
    x_int = zeros(size(Resp_r,1),npoints);
    
    for i=1:size(Resp_r,1)
    
        %estimate per output
        LB = min(Resp_r(i,1)-5*Resp_r(i,2),Resp_s(i,1)-5*Resp_s(i,2)); %lower bound for integral
        UB = max(Resp_r(i,1)+5*Resp_r(i,2),Resp_s(i,1)+5*Resp_s(i,2)); %upper bound for integral 
        x_int = linspace(LB,UB,npoints); % integration points
        
        pdf_r = normpdf(x_int,Resp_r(i,1),Resp_r(i,2));
        pdf_s = normpdf(x_int,Resp_s(i,1),Resp_s(i,2));
        
        H_s(i) = 1-trapz(x_int,sqrt(pdf_s.*pdf_r)); % trapezoidal integration
        
    end
    
    % average squared Hellinger distance per output
    H_s = mean(H_s); 


end

