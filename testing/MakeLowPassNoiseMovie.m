%% MakeLowPassNoiseMovie.m
% written by G, M. Boynton, 2017, at the University of Washington 
% 
% This script creates and displays a "movie" of random noise that has been
% filtered so that the power spectrum of each spatial frequency *and* temporal
% frequency component follows the function 1/f^n (where f is the
% frequency).
% 

%% Parameters
powerx = -1;  %power of 1/f^n noise for space
powert = -1;  %power of 1/f^n noise for time
nx = 201;    %square image size (pixels)
nt = 200;   % number of frames


%% initialize movie frames, starting with Gaussian noise
img = zeros(nx,nx,nt);
noise = randn(nx,nx,nt);


%% Filter each frame spatially: 
for i=1:nt
    image0 = noise(:,:,i);
    image0 = image0-mean(image0(:));        % de-mean
    F = myfft2(image0);                     % get 2D spectrum
    mypower = F.sf.^powerx;                 % make coefficents proprtional to f^n
    F.amp = F.amp.*mypower;                 % modify amplitudes
    img(:,:,i) = myifft2(F);
end

%% Filter each pixel temporally: 
for i=1:nx
    for j=1:nx
        F = complex2real(fft(squeeze(img(i,j,:))));
        F.amp = F.amp.*F.freq.^powert;
        img(i,j,:) = real(ifft(real2complex(F)));
    end
end


%% scale the image for proper display:
mx = max(abs(img(:)));
img = 256*(img/mx+1)/2;

%% show the movie

frameRate = 30; % in Hz, frames per second 
nReps = 2; %number of times to loop the movie

figure(1)
clf
h = image(img(:,:,1));
colormap(gray(256))
axis equal
axis off
hold on
set(gca,'Position',[0,0,1,1]);

for j=1:nReps
    tic
    for i=1:nt
        set(h,'CData',img(:,:,i));
        while toc<i/frameRate
        end
        drawnow
    end
end





