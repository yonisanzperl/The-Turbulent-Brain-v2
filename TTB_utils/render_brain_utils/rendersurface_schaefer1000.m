function rendersurface_schaefer1000(schaefer1000vector,rangemin,rangemax, inv,clmap,surfacetype, fig_name)
% script for rendering a schaefer1000vector
%   ML Kringelbach April 2020
% schaefer1000vector : the schaefer1000vector values to be rendered
%
% rangemin, rangemax : limits for colorscheme
%
% inv:
%  0 colormap, interp
%  1 flip colormap, interp
%  2 colormap, only three colours
%
% clmap : from othercolor, default 'Bu_10'
%
% surfacetype: 
%  1 midthickness
%  2 inflated (default)
%  3 very inflated

% colormaps and tight subplots

addpath(genpath([pwd '/TTB_utils/render_brain_utils/']))


if ~exist('rangemin','var')
     rangemin=min(schaefer1000vector);
end

if ~exist('rangemax','var')
     rangemax=max(schaefer1000vector);
end

if ~exist('inv','var')
     inv=max(schaefer1000vector);
end

if ~exist('clmap','var')
    clmap='Bu_10';
end

if ~exist('surfacetype','var')
     surfacetype=2; % default is inflated
end


% make space tight
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

% load the different views
base=[pwd '/TTB_utils/render_brain_utils/'];;
display_surf_left=gifti([base 'ParcellationPilot.L.inflated.32k_fs_LR.surf.gii']);
display_surf_right=gifti([base 'ParcellationPilot.R.inflated.32k_fs_LR.surf.gii']);

basedir=[pwd '/TTB_utils/render_brain_utils/'];;
glassers_L=gifti([basedir 'Glasser360.L.mid.32k_fs_LR.surf.gii']);
glassersi_L=gifti([basedir 'Glasser360.L.inflated.32k_fs_LR.surf.gii']);
glassersvi_L=gifti([basedir 'Glasser360.L.very_inflated.32k_fs_LR.surf.gii']);
glassersf_L=gifti([basedir 'Glasser360.L.flat.32k_fs_LR.surf.gii']);
glassers_R=gifti([basedir 'Glasser360.R.mid.32k_fs_LR.surf.gii']);
glassersi_R=gifti([basedir 'Glasser360.R.inflated.32k_fs_LR.surf.gii']);
glassersvi_R=gifti([basedir 'Glasser360.R.very_inflated.32k_fs_LR.surf.gii']);
glassersf_R=gifti([basedir 'Glasser360.R.flat.32k_fs_LR.surf.gii']);

switch surfacetype
    case 1
        display_surf_left=glassers_L;
        display_surf_right=glassers_R;
    case 2
        display_surf_left=glassersi_L;
        display_surf_right=glassersi_R;
    case 3
        display_surf_left=glassersvi_L;
        display_surf_right=glassersvi_R;
end;

sl = display_surf_left;
sr = display_surf_right;


%base='/Users/mortenk/Documents/MATLAB/schaefer/';
base = [pwd '/TTB_utils/render_brain_utils/'];
atlas_l=gifti([base 'Schaefer1000_L.func.gii']);
atlas_r=gifti([base 'Schaefer1000_R.func.gii']);
label=ft_read_cifti([base 'Schaefer2018_1000Parcels_7Networks_order.dscalar.nii']);

vl=atlas_l;
vr=atlas_r;


% % replace the left hemisphere labels with values from dkt(1:34)
% % dk46_L.labels contains 34 elements
% left hemisphere
for i=1:500
    idx=find(atlas_l.cdata==i);
    vl.cdata(idx)=schaefer1000vector(i);
end;

% right hemisphere
for i=501:1000
    idx=find(atlas_r.cdata==i);
    vr.cdata(idx)=schaefer1000vector(i);
end;

% remove -1 from labels
% idx{i}=find(label_L.cdata==-1);
% vl.cdata(idx{i})=0;
% idx{i}=find(label_R.cdata==-1);
% vr.cdata(idx{i})=0;


%% rendering

    % create figure
    hfig = figure;
    hfig.Name = fig_name;
    set(gcf, 'Position',  [40, 40, 400, 600]);
    
    subplot(3,2,1); %left hemisphere side view
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(1) = patch(ax2,'Faces',sl.faces,'vertices',sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(-90,0);
    camlight;
    lighting gouraud;
    material dull;


    subplot(3,2,3); %left hemisphere midline
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(1) = patch(ax2,'Faces',sl.faces,'vertices',sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(90,0)
    camlight;
    lighting gouraud;
    material dull;

    subplot(3,2,4); %right hemisphere side view
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sr.faces,'vertices',sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(-90,0)
    camlight;
    lighting gouraud;
    material dull;
 
    subplot(3,2,2); %right hemisphere midline
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sr.faces,'vertices',sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(90,0)
    camlight;
    lighting gouraud;
    material dull;

    % flatmaps
    sl=glassersf_L;
    sr=glassersf_R;
    
    subplot(3,2,5); %left hemisphere flat
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sl.faces,'vertices',sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(0,90)
    camlight;
    lighting gouraud;
    material dull;
    %colorbar('southoutside');

    subplot(3,2,6); %right hemisphere flat
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sr.faces,'vertices',sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(0,90)
    camlight;
    lighting gouraud;
    material dull;
    %colorbar('southoutside');
    colorbar('Location', 'southoutside', 'Position', [.25 .05 .5 .05])
    switch inv
        case 0
            % use selected colormap (interpolated to 64 values)
            c=othercolor(clmap);
            colormap(c)            
        case 1
            % flip colormap (interpolated to 64 values)
            c=flipud(othercolor(clmap));
        case 2
            % use specialised version with only three values
            c=othercolor(clmap,3);
            % this is for the second value
            c(2,1)=0.70;c(2,2)=0.70;c(2,3)=0.70; 
        case 3
            % use specialised version with only four values
            c=othercolor(clmap,4);
            c(2,1)=0.70;c(2,2)=0.70;c(2,3)=0.70; 
    end;
    % neutral brain colour: grey
    c(1,1)=0.98;c(1,2)=0.98;c(1,3)=0.98; 
    colormap(c)
    %colorbar

