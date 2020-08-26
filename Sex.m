function [W] = Sex(mu,rsex,rrsex,fr)
   
    M=20; % M is the number of discrete mitochondria per cell
    n=1;
    x=4;  % x is the exponent in the fitness function. x>1 makes it concave which we want.
    s=1;  % selection strength

    P=zeros(M+1,n);
    p = rand([M+1 1]);
    P(:,1)=1*p./sum(p);

    PR=zeros(M+1,n);
   
    % Life cycle operators
    % Mitochondrial mutation
    U=zeros(M+1,M+1);
    for i=0:M
        for j=0:M
            U(i+1,j+1)=sum(binopdf(i-j,M-j,mu));
        end
    end
   
    w=zeros(M+1,n);
    for i=0:M
        w(i+1,1)=1 - s*( i / M )^x;
    end
    
    % Resampling at Meiosis 1
    L1=zeros(2*M+1,2*M+1);
    for i=0:2*M
        for j=0:2*M
            L1(i+1,j+1)=hygepdf(i,4*M,2*j,2*M);
        end
    end
    
    % Resampling at Meiosis 2
    L2=zeros(M+1,2*M+1);
    for i=0:M
        for j=0:2*M
            L2(i+1,j+1)=hygepdf(i,2*M,j,M);
        end
    end

    AD=zeros(0,1);
    for i=0:2*M
        for j=0:M
            AD(i+1,j+1)=binopdf(i,2*M,j/M);
        end
    end
    % Now iterate2the life cycle for many many generations k
    for k=1:2000
         if k==1
             PR=fr*P;
             P=P-fr*P;
         end
        
        % Mutate mitochondria from a to A
        
        P = U*P;
        PR = U*PR;
        % Selection
        PPR = ([w,w].*[P,PR]) / sum([w,w].*[P,PR],'All');
        P = PPR(:,1);
        PR = PPR(:,2);
        
        % Sexual distributions S
        % Z - P+PR, PR+PR
       
        Z(:,1)=conv(P,P);
        Z(:,2)=conv(P,PR)*2;
        Z(:,3)=conv(PR,PR);
                
        Zsum = sum(Z,'All');
        if Zsum==0
            Zsum=1;
        end
        Z = Z/Zsum;
        
        
        S=L2*L1*( Z(:,1)+0.5*Z(:,2));
        SR=L2*L1*( 0.5*Z(:,2)+Z(:,3) );
                      

        % asexual and sexual combined
        P=L2*L1*AD*P*(1-rsex)+S*rsex;
        PR=L2*L1*AD*PR*(1-rrsex)+SR*rrsex;
                        
        % if we need to keep mutant frequency constant to measure fitness
        % advantage, set it to fr here:
        P = P/sum(P,'All')*(1-fr);
        PR = PR/sum(PR,'All')*(fr);
        
        % record mean fitness
        %W(k) = sum(w.*P,'All')/sum(P,'All'); 
        W(k) = sum([w,w].*[P,PR],'All');
        WR(k)= sum(w.*PR,'All')/sum(PR,'All');
       
        %W(k)=sum(PR,'All');
    end
   
    W = WR(end)-W(end);
    
    Pmito = [P/sum(P),PR/sum(PR)];
   
    

end
    