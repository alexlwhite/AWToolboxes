%must first have loaded in vtc data with a command like: 
%vtc = BVQXfile(vtcFile);
% 
slice = 40; %pick a slice to view
figure(1); set(gcf,'units','pixels','pos',[5 5 800 700]);
nt = size(vtc.VTCData,1);
for t = 1:nt
    figure(1); clf;
    img = squeeze(vtc.VTCData(t,:,slice,:));
    imagesc(img); colormap gray;
    title(num2str(t));
    pause(0.05);
end
