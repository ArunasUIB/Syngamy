function [W] = SexUnequalDivision(M,rsex,rrsex,fr)
    n=1;
    x=2;
    P=zeros(M+1,n);
    p = rand([M+1 1]);
    P(:,1)=1*p./sum(p);

    PR=zeros(M+1,n);
   
%     w=ones(M+1,n);
%     w(1:2,1)=0.1;
%     w(M+1,1)=0;
    w=zeros(M+1,1);
    for i=0:M
        w(i+1,1)=1 - (1- i/M)^x;
    end
         
    AD=zeros(0,1);
    for i=0:2*M
        for j=0:M
            AD(i+1,j+1)=(i==2*j);
        end
    end
       
    B=zeros(0,1);
    for i=0:M
        for j=0:2*M
            B(i+1,j+1)=binopdf(i,j,1/2);
        end
    end        
    
    % Now iterate the life cycle for many many generations k
    for k=1:1000
         if k==1
             PR=fr*P;
             P=P-fr*P;
         end
        
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
      
        S=B*( Z(:,1)+0.5*Z(:,2));
        SR=B*( 0.5*Z(:,2)+Z(:,3) );
                      

        % asexual and sexual combined
        P=B*AD*P*(1-rsex)+S*rsex;
        PR=B*AD*PR*(1-rrsex)+SR*rrsex;
        
        
                        
        % if we need to keep mutant frequency constant to measure fitness
        % advantage, set it to fr here:
        P = P/sum(P,'All')*(1-fr);
        PR = PR/sum(PR,'All')*(fr);
        
        % record mean fitness
        %W(k) = sum(w.*P,'All')/sum(P,'All'); 
        W(k) = sum([w,w].*[P,PR],'All');
        WR(k)= sum(w.*PR,'All')/sum(PR,'All');
        Wmut(k)=sum(PR,'All');
        
    end
    
    W = WR(end)-W(end);
    
    Pmito = [P/sum(P),PR/sum(PR)];
    mt=0:M;
    sum(P.*mt.')/sum(P,'All');
    sum(PR.*mt.')/sum(PR,'All');
    
    plot(0:M,Pmito)
   
    

end
    