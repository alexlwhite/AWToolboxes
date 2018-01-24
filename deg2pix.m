function pix=deg2pix(MO_WIDE,MO_PHYS,ABSTAND,DEGS)

phi = atan2(1,ABSTAND)*180/pi;  % wie viel Grad entsprechen 1 cm auf dem Bildschirm?
help = MO_WIDE/MO_PHYS;         % MO_WIDE Pixel (x-Achse) entsprechen 35 cm => help gibt Verhaeltnis von Pixeln und cm an
phi = phi / help;               % DEG_PER_PIX
pix = DEGS/phi;