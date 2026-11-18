function write_frame(fig,k,kmin,gifname,delayT)
% Writes a frame of an animation to a .gif

frame = getframe(fig); % Capture the frame
img = frame2im(frame);
[img,cmap] = rgb2ind(img,256);
if k == kmin
    imwrite(img,cmap,gifname,'gif','LoopCount',Inf,'DelayTime',delayT);
else
    imwrite(img,cmap,gifname,'gif','WriteMode','append','DelayTime',delayT);
end

end