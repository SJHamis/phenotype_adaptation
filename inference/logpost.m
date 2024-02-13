function lp = logpost(expars,model,Ydata,lower,upper)
    
    % Extended support
    exlower = [0,lower];
    exupper = [10,upper];

    % Prior support function
    insupport = @(expars) all((expars > exlower) & (expars < exupper));


    if insupport(expars)
        lp = loglike(expars,model,Ydata);
    else
        lp = -Inf;
    end

end

    % Likelihood function
    function ll = loglike(expars,model,Ydata)
    
        % Noise std
        std = expars(1);
    
        % Model parameters
        pars = expars(2:end);
        
        % Model prediction
        Ymdl = model(pars);
    
        % Log likelihood
        eps = log(Ymdl) - log(Ydata);
        ll = sum(0.5 * (eps(:) / std).^2 + 1 / (std * sqrt(2 * pi)));
    
    end
