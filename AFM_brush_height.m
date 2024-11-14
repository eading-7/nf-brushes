import IBWread.*
import readIBWheaders.*
import readIBWbinheader.*
%Read raw data from .ibw filetype with:
%Jakub Bialek (2024). Igor Pro file format (ibw) to matlab variable (https://www.mathworks.com/matlabcentral/fileexchange/42679-igor-pro-file-format-ibw-to-matlab-variable),
%MATLAB Central File Exchange. Retrieved November 13, 2024.

%This file extracts protein brush height from raw .ibw force curve files, exported from Oxford Instruments AsylumResearch AFM Software.
%Code written by Ding, E. A. UC Berkeley, Department of Chemical and Biomolecular Engineering.

finalresults = [];
address = ''; %data path
foldernames = {''}; %names of folders containing raw .ibw files

for foldernum = 1:size(foldernames,2)
    clearvars -except foldernames foldernum finalresults address
    
    %Import data
    foldername = char(foldernames(foldernum));
    d = dir([address,foldername,'\*.ibw']);
    nfiles = length(d);
    alldata(nfiles) = struct();
    for jj = 1:nfiles
        alldata(jj).data = IBWread([d(jj).folder,'\',d(jj).name]).y;
    end
    
    %Read the cantilever spring constant for the force map, in N/m, from the metadata
    metadata = splitlines(IBWread([d(1).folder,'\',d(1).name]).WaveNotes);
    metadata = cellfun(@(x) regexp(x,':','split','once'),metadata(1:end - 1),'UniformOutput',false);
    metadata = vertcat(metadata{:});
    k_c = metadata{strcmp(metadata,'SpringConstant'),2};
    k_c = str2double(k_c); %N/m
    
    n_discard = 0;
    height = nan(nfiles,1);
    
    close('all');
    disp(foldername);
    figure(1);
    ylabel('Force (N)');
    xlabel('Separation (m)');
    title('Accepted Curves');
    hold on
    figure(2);
    ylabel('Force (N)');
    xlabel('Z position (m)');
    title('Discarded Curves (force distance)');
    hold on
    
    for I = 1:nfiles
        %Only use data from the approach curve
        %column 1: Z position in ramp, m
        %column 2: cantilever deflection, m
        [maxpoint,approach] = max(alldata(I).data(:,1));
        a = alldata(I).data(1:approach,:);
        
        nr = length(a); nbase = floor(0.01*nr)+1;
        nend = 25;      nsep = floor(0.25*nr)+1; %fit region sizes
        tolcon = 3; tolsep = 3; %fit tolerances
        
        asep = a(nbase:nbase+nsep,:);
        
        [bc,bc_s] = polyfit(asep(:,1),asep(:,2),1);
        a(:,2) = (a(:,2) - bc(1)*a(:,1)-bc(2))*k_c; %baseline-correct the deflection curve and convert to force, N
        a(:,1) = a(:,1) - a(:,2)/k_c; %get ramp position of tip, m
        
        asep = a(nbase:nbase+nsep,:);
        sigsep = norm(asep(:,2))/sqrt(nsep);
        
        %Reject curves with surface contamination
        smooth = 0;
        if sigsep<0.01*max(a(:,2))
            dxdy = diff(movmean(a(:,1),25))./diff(movmean(a(:,2),25));
            acontact = a(nr-nend-1:nr-1,:);
            [p,s] = polyfit(acontact(:,2),acontact(:,1),1); %flip fit because near vertical: x = p(1)y + p(2)
            sigcon = s.normr/sqrt(s.df);
            dxdypos = (dxdy > -1);
            islope = floor(0.5*nr);
            while any(dxdypos(islope:islope+12)) && islope < nr-13
                islope = islope + 1;
            end
            if abs(a(islope,1)-p(1)*a(islope,2)-p(2)) < sigcon*10 || range(a(islope:islope+25,2)) < 1e-10
                smooth = 1;
            end
        end
        
        if smooth
            %%Fit contact region
            %Calculate deviation from linear fit
            contactfit = polyval(p,a(:,2));
            devcon = a(:,1)-contactfit;
            
            %Approach the contact region until consistently within tolerance of the linear contact region fit
            icon = floor(nr*0.5);
            while any(abs(devcon(icon:icon+2))>tolcon*sigcon)
                icon = icon + 1;
            end
            z_contact = a(icon,1);
            
            %%Fit baseline separation
            %Approach the baseline until consistently within tolerance of the baseline
            isep = nr - nend;
            while any(a(isep-2:isep,2)>tolsep*sigsep)
                isep = isep - 1;
            end
            z_sep = a(isep,1);
            
            height(I) = z_contact - z_sep; %brush height, m
            
            %%Plot all curves for visual verification of data quality and processing
            %Align contact points for plotting
            z_aligned = z_contact-a(:,1);
            
            figure(1);
            plot(z_aligned, a(:,2),'b'); %data
            plot(z_contact-contactfit(nr-nend:nr),a(nr-nend:nr,2),'r','LineWidth',2); %contact fit
            plot(z_contact-z_sep,a(isep,2),'ko','MarkerSize',10); %separation point
            plot(0,a(icon,2),'ro','MarkerSize',10); %contact point
            if z_contact - z_sep < 0 
                plot(z_contact-z_sep,a(isep,2),'gx','MarkerSize',10);
            end
        else
            %For rejected curves, plot force vs Z position
            figure(2);
            plot(a(:,1),a(:,2));
            n_discard = n_discard+1;
            disp(I);
        end
    end
    
    %remove the placeholder values for any rejected curves
    height = rmmissing(height);
    %report on data quality
    disp('Number of Used Ramps:');
    disp(nfiles-n_discard);
    disp('Number of Unused Ramps:');
    disp(n_discard);
    
    saveas(figure(1),[address,foldername,'.png'])
    saveas(figure(2),[address,foldername,'rejected.png'])
    
    finalresults{foldernum} = height;
end

writecell(finalresults','finalresults.csv');