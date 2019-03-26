%MakeLowPassColorNoiseMovie.m
%Adapted by ALW from GMBs makeLowPassNoiseMovie, to have independent noise
%in the 3 color channels. 

powerx = -1;  %power of 1/f^n noise for space
powert = -3;  %power of 1/f^n noise for time
nx = 601;    %square image size (pixels)
nt = 200;   % number of frames
nc = 3;      %number of color channels

img = zeros(nx,nx,nc,nt);
noise = randn(nx,nx,nc,nt);


%Generate a sequence of 1/f^n images in space
for c=1:nc
    for i=1:nt
        image0 = noise(:,:,c,i);
        image0 = image0-mean(image0(:));
        F = myfft2(image0);                     % get 2D spectrum
        mypower = F.sf.^powerx;                 % make coefficents proprtional to f^n
        F.amp = F.amp.*mypower;                 % modify amplitudes
        img(:,:,c,i) = myifft2(F);
    end
end

% Now do the same in time for every pixel
for i=1:nx
    for j=1:nx
        for c=1:nc
            F = complex2real(fft(squeeze(img(i,j,c,:))));
            F.amp = F.amp.*F.freq.^powert;
            img(i,j,c,:) = real(ifft(real2complex(F)));
        end
    end
end


% scale the image
mx = max(abs(img(:)));
img = (img/mx+1)/2; %for color, don't multiply by 256

%%
% show the movie

frameRate = 30;
nReps =3 ;

figure(1)
clf
h = image(squeeze(img(:,:,:,1)));
%colormap(gray(256))
axis equal
axis off
hold on
set(gca,'Position',[0,0,1,1]);

for j=1:nReps
    tic
    for i=1:nt
        set(h,'CData',squeeze(img(:,:,:,i)));
        while toc<i/frameRate;
        end
        drawnow
    end
end





