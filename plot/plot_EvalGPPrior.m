function plot_EvalGPPrior(t, y, covy_emp, covy_hp)
% compares empirical (normal) distribution to modeled normal distribution 

% INPUT:
% t (1-by-nsteps) = time
% y (ny-by-nsteps) = measurements
% covy_emp (ny-by-ny) = empirical covariance of measured states
% covy_hp (ny-by-ny) = modeled covariance of measured states

% OUTPUT:
% figure

    
    % extract nonzero rows
    idx = find(all(y(:,2:end),2)); 
    y_nz = y(idx,:);
    blue = [86 180 233]/255;
    orange = [230 159 0]/255;
        
    ny = size(y_nz,1);
    nt = size(y,2);
    std_emp = sqrt(diag(covy_emp));
    std_hp = sqrt(diag(covy_hp));
    
    
    figure; set(gcf, 'Position', [100, 100, 700, 200*ny]);
    for i = 1:ny
    
        % time series
        sp1 = subplot_tight(ny,2,2*i-1,[0.07 0.08]);
        sp1.Position = sp1.Position + [0 0 0.2 0];
            
            % plot measurement series
            plot(t, y_nz(i,:), 'Color', [0.2 0.2 0.2]); hold on
        
            % plot empirical distribution
            ytop = ones(1,nt)*3*std_emp(i); ybot = -ytop;
            plot(t, ytop, 'Color', blue, 'Linewidth', 2);
            plot(t, ybot, 'Color', blue, 'Linewidth', 2);
        
            % plot modeled prior
            ytop = ones(1,nt)*3*std_hp(i); ybot = -ytop;
            plot(t, ytop, '--', 'Color', orange, 'Linewidth', 2);
            plot(t, ybot, '--', 'Color', orange, 'Linewidth', 2);
            
            % labels
            ylabel(sprintf("channel %g", i));
            if i == ny
                xlabel('$t$', 'interpreter', 'latex');
            elseif i == 1
                legend("obs.", "$\pm3\sigma$  (obs.)", "", "$\pm3\hat{\sigma}_0$ (prior)", 'interpreter', 'latex', 'Location', 'southwest');
            end
        
        
        % pdfs
        sp2 = subplot_tight(ny,2,2*i,[0.07 0.08]);
        sp2.Position = sp2.Position + [0.2 0 -0.2 0];

            plotrng = min(y_nz(i,:)):max(y_nz(i,:))/200:max(y_nz(i,:));
            pd = fitdist(y_nz(i,:)', 'Normal');
            pdf_emp = pdf(pd, plotrng);
            pdf_hp = normpdf(plotrng,0,std_hp(i));
            [counts, centers] = histcounts(y_nz(i,:), 100, 'Normalization', 'pdf');
            
            barh(centers(1:end-1), counts, 'FaceColor', [0.2 0.2 0.2]); hold on
            plot(pdf_emp, plotrng, 'Color', blue, 'Linewidth', 2); 
            plot(pdf_hp, plotrng, '--', 'Color', orange, 'Linewidth', 2);
            if i == ny; xlabel('$PDF [-]$', 'interpreter', 'latex'); end
    

    end


end