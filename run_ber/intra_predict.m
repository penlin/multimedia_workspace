function mse = intra_predict(beta,file,gamma)

im = double(rgb2gray(imread(file)));
[h w ] = size(im);
img_bp = img2bp_frame_mex(im);
PXL = 8;
order = 4.^(-[-7:0]);
x = csvread('xx.txt');
y = csvread('yy.txt');

fr = interp1(x,y,gamma,'cubic')
%fr_matrix = zeros(h,w,PXL);
frr = zeros(1,8);

for i = 1:h,
    for j = 1:w,
        if i == 1 || i== h || j == 1 || j == w,
           %fr_matrix(i,j,:) = fr(:);
           frr = frr + fr;
           continue; 
        end
%        for t_lvl = 1:PXL,
        n = squeeze(img_bp(i+1,j,:)+img_bp(i-1,j,:)+img_bp(i,j+1,:)+img_bp(i,j-1,:))';
        %fr./(fr+(1-fr).*exp((2*squeeze(img_bp(i,j,:))'-1).*beta.*(4-2*n)))
        frr = frr + fr./(fr+(1-fr).*exp((2*squeeze(img_bp(i,j,:))'-1).*beta.*(2*n-4)));
%        end
    end
end

frr = frr/h/w

mse = order*frr';

end