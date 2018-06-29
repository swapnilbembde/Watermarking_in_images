clc;
clear all;
warning off;
chos=0;
possibility=10;



while chos~=possibility,
    chos=menu('Digital image watermarking','Select Cover Image','Select watermarking image','Show Spatial Domain watermarking','Show Spectral Watermarking using DCT','Show Spectral Watermarking using DFT','show 3-level coverimage','show 3-level watermarkimage','show watermarked image','show extracted image','exit');
    if chos==1
        [fname pname]=uigetfile('*.jpg','select the Cover Image');
        %eval('imageinput=imread(fname)');
        imageinput=imread(fname);
        

A=rgb2gray(imageinput);
figure
imshow(A,'DisplayRange',[]), title('Selected cover image');
    end
    if chos==2
    [fname pname]=uigetfile('*.jpg','select the Watermark');
    imageinput=imread(fname);
    W2=rgb2gray(imageinput);
    W = imresize(W2,size(A));
figure
imshow(W,'DisplayRange',[]), title('Selected Watermarking image');
    end
    if chos==3
        [B1,B2,B3,B4,B5,B6,B7,B8] = bitplane(A);
        [D1,D2,D3,D4,D5,D6,D7,D8] = bitplane(W);
        Final = B1 * 2^7 + B2 * 2^6 + B3 * 2^5 + B4 * 2^4 + D1 * 2^3 + D2 * 2^2 + D3 * 2^1 + D4 * 2^0;
        figure
        imshow(Final),title('Watermarked Image');
    end
    if chos==4
            img_gray_dct = dct2(A);
            a = size(img_gray_dct, 1);
            b = size(img_gray_dct, 2);

            dct_a = reshape(img_gray_dct, 1, a*b);

            water_gray = W;
            water_small = imresize(water_gray, [floor(3*a/4), floor(3*b/4)]);
            [max_val, ind1] = max(water_small, [], 1);
            [max_water, ind2] = max(max_val, [], 2);
            water_small = double(water_small)/double(max_water);
            disp(size(water_small));
            figure;
            imshow(uint8(water_small*double(max_water)));

            water_a = size(water_small, 1);
            water_b = size(water_small, 2);

            water_dct = dct2(water_small);
            water_dct1 = reshape(water_dct, 1, water_a*water_b);

            k = water_a*water_b;

            [real_dct_sort, index] = sort(abs(dct_a), 'descend');
            disp(size(real_dct_sort));

            dct_max = [];
            dct_max = real_dct_sort(1, 2:k);
            disp(size(dct_max));

            alpha = 0.0001;

            watermarked = [];
            dct_c = dct_a;
            for j=1:k-1
                dct_c(index(j)) = dct_a(index(j)) + (alpha*abs(dct_a(index(j)))*water_dct1(j));
            %     watermarked(j) = dft_max(j) + (alpha*abs(dft_max)*water(j));
            end

            dct_b = reshape(dct_c, a, b);
            watermarked_img = idct2(dct_b);
            figure;
            imshow(uint8(watermarked_img));

            imwrite(uint8(watermarked_img), 'watermarked_image.png');
            img1 = imread('watermarked_image', 'png');
            dct_b = dct2(img1);

            for j=1:k-1
                extracted(j) = (dct_b(index(j)) - dct_a(index(j)))/(alpha*abs(dct_a(index(j))));
            end

            for c = 1:10^2
                rand_seq = randn(1,k-1);
                if(c == 25)
                    rand_seq = water_dct1(1, 1:k-1);
                end
                SIM_p(c) = sum(extracted.*rand_seq);
end

plot(abs(SIM_p));
    end
    if chos==5
        img_gray_dft = fftshift(fft2(A));
        a = size(img_gray_dft, 1);
        b = size(img_gray_dft, 2);

        dft_a = reshape(img_gray_dft, 1, a*b);

        water_img = W;
        water_gray = W;
        water_small = imresize(water_gray, [floor(3*a/4), floor(3*b/4)]);
        [max_val, ind1] = max(water_small, [], 1);
        [max_water, ind2] = max(max_val, [], 2);
        water_small = double(water_small)/double(max_water);
        disp(size(water_small));
        figure;
        imshow(uint8(water_small*double(max_water)));

        water_a = size(water_small, 1);
        water_b = size(water_small, 2);

        water_dft = fftshift(fft2(water_small));
        water_dft1 = reshape(water_dft, 1, water_a*water_b);

        k = water_a*water_b;

        [real_dft_sort, index] = sort(abs(dft_a), 'descend');
        disp(size(real_dft_sort));

        dft_max = [];
        dft_max = real_dft_sort(1, 2:k);
        disp(size(dft_max));

        alpha = 0.001;

        watermarked = [];
        dft_c = dft_a;
        for j=1:k-1
            dft_c(index(j)) = dft_a(index(j)) + (alpha*abs(dft_a(index(j)))*water_dft1(j));
        end

        dft_b = reshape(dft_c, a, b);
        watermarked_img = ifft2(ifftshift(dft_b));
        figure;
        imshow(uint8(watermarked_img));

        imwrite(uint8(watermarked_img), 'watermarked_image.jpg', 'Mode', 'lossy', 'Quality', 75);
        img1 = imread('watermarked_image', 'jpg');
        figure;
        imshow(img1), title('Watermarked after compression');
        dft_b = fftshift(fft2(img1));

        for j=1:k-1
            kaala_paani(j) = (dft_b(index(j)) - dft_a(index(j)))/(alpha*abs(dft_a(index(j))));
        end

        for c = 1:10^2
            rand_seq = randn(1,k-1);
            if(c == 25)
                rand_seq = water_dft1(1, 1:k-1);
            end
            SIM_p(c) = sum(kaala_paani.*rand_seq);
        end

        plot(abs(SIM_p));
    end
    if chos==6
        P1=im2double(A);
P=imresize(P1,[2048 2048]);
    [F1,F2]= wfilters('haar', 'd');
[LL,LH,HL,HH] = dwt2(P,'haar','d');
[LL1,LH1,HL1,HH1] = dwt2(LL,'haar','d');
[LL2,LH2,HL2,HH2] = dwt2(LL1,'haar','d');
%figure(2)
imshow(LL,'DisplayRange',[]), title('3 level dwt of cover image');
%title('extracted watermark')
    end
    if chos==7
        watermark=im2double(W2);
watermark=imresize(watermark,[2048 2048]);
%figure(3)
[WF1,WF2]= wfilters('haar', 'd');
[L_L,L_H,H_L,H_H] = dwt2(watermark,'haar','d');
[L_L1,L_H1,H_L1,H_H1] = dwt2(L_L,'haar','d');
[L_L2,L_H2,H_L2,H_H2] = dwt2(L_L1,'haar','d');
imshow(L_L2,'DisplayRange',[]), title('3 level dwt of watermark image');
    end
    if chos==8
        k_h=1; q_h =0.009;
        Watermarkedimage=k_h*LL2+q_h*L_L2;
        %computing level-1 idwt2
Watermarkedimage_level1= idwt2(Watermarkedimage,LH2,HL2,HH2,'haar');
%figure(5)

%computing level-2 idwt2
Watermarkedimage_level2=idwt2(Watermarkedimage_level1,LH1,HL1,HH1,'haar');
%figure(6)

%computing level-3 idwt2
Watermarkedimage_final=idwt2(Watermarkedimage_level2,LH,HL,HH,'haar');
%figure(7)
imshow(Watermarkedimage_final,'DisplayRange',[]), title('Watermarkedimage final')
    end
     if chos==9
        [F11,F22]= wfilters('haar', 'd');
[a b c d]=dwt2(Watermarkedimage_final,'haar','d');
[aa bb cc dd]=dwt2(a,'haar','d');
[aaa bbb ccc ddd]=dwt2(aa,'haar','d');

recovered_image=aaa-k_h*LL2;
%figure(8)
imshow(recovered_image,[]);
%title('extracted watermark')
    end
end
