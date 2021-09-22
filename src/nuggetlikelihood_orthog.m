function[val] = nuggetlikelihood_orthog(x,Phi_S_Phi,trS,n)

        l = size(Phi_S_Phi,1);
        
        Phi_S_Phi_over_nugsq = Phi_S_Phi/(x(1)^2);

        for j=1:l
            Phi_S_Phi_over_nugsq(:,j) = Phi_S_Phi_over_nugsq(:,j) * 1/(x(2)+1/(x(1)));
        end

        
        logdetpart = l*log(x(2)+1/(x(1))) - l*log(x(2));
        tracepart  = sum(diag( Phi_S_Phi_over_nugsq  ));

        constantpart= n*log(x(1)) + trS/(x(1));

        
        val =  logdetpart - tracepart + constantpart;
end

