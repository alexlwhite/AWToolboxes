%testing HSV colors 
nhues=100;
nsats=50;

huesL=linspace(0,1,nhues); 
satsL=linspace(0,1,nsats); 

val=.85; 

hues=repmat(huesL,nsats,1); 
sats=repmat(satsL',1,nhues); 
vals=ones(size(hues))*val; 

rgbsL=hsv2rgb([hues(:) sats(:) vals(:)]); 

%reshape 
rgbs=zeros(nsats,nhues,3);
for gi=1:3
    rgbs(:,:,gi)=reshape(rgbsL(:,gi),nsats,nhues);
end

figure; 
%imshow(rgbs)

figure; hold on;
for si=1:nsats
    for hi=1:nhues
        plot(satsL(si), huesL(hi), '.','Color',squeeze(rgbs(si,hi,:)),'MarkerSize',30);
    end
end