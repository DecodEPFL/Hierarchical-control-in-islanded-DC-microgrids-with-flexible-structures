function [Ac, Buc, Mdc]=MPC_matmaker(A,Bu,Md,N)


[ra,ca]=size(A);
Ac=A;


    for i=2:N  
        
        Ac=[Ac;A^(i)];
    
    end
    Ac=[eye(ra);Ac];
    

[rb,cb]=size(Bu);

Buc=[];

    for i=1:N
        
        Bucr=[];
     
        for j=1:N
            
            if(i>=j)
                Bucr=[Bucr,(A^(i-j))*Bu];
            else
                Bucr=[Bucr,zeros(ra,cb)];
            end
        end
        
        Buc=[Buc; Bucr];
    end
    Buc=[zeros(rb,size(Buc,2)); Buc];
    
    

[rbw,cbw]=size(Md);
Mdc=[];

    for i=1:N
        
        Mdcr=[];
     
        for j=1:N
            
            if(i>=j)
                Mdcr=[Mdcr,(A^(i-j))*Md];
            else
                Mdcr=[Mdcr,zeros(ra,cbw)];
            end
        end
        
        Mdc=[Mdc; Mdcr];
    end 
    Mdc=[zeros(rbw,size(Mdc,2)); Mdc];

    
end
    

