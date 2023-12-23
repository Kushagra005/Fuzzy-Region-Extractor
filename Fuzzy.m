clear all;
a = imread('6.jpg');
tic
dim = 100; %query image size
diml = 2100; %database image size
maxtr = [];
maxtg = [];
maxtb = [];

%% Calculating the maximum count for each channel of all windows of the scene image
for il = 1:dim:diml
    for i2 = 1:dim:diml
    im = imcrop(a,[i2 il dim-1 dim-1]);
    imr = im(:,:,1);
    img = im(:,:,2);
    imb = im(:,:,3);
    %histogram for red channel%
    [count,xlr]=imhist(imr);
    countimr=count;
    %histogram for green channel%
    [count,xlg]=imhist(img);
    countimg=count;
    %histogram for blue channel%
    [count,xlb]=imhist(imb);
    countimb=count;
    %maximum count of red channel in each window of the scene image%
    maxllr=max(countimr);
    maxtlr = [maxtr,maxllr];
    maxtr = maxtlr;
    %maximum count of Green channel in each window of the scene image%
    maxllg=max(countimg);
    maxtlg=[maxtg,maxllg];
    maxtg=maxtlg;
    %maximum count of Blue channel in each window of the scene image%
    maxllb = max(countimb);
    maxtlb=[maxtb,maxllb];
    maxtb=maxtlb;
    end
end

%maximum count among all the windows in each channel(normalizing factor)%
mx2r = max(maxtlr);
mx2g = max(maxtlg);
mx2b = max(maxtlb);

%membership value calculation for the query image, c(query image) are
%calculated by the normalizing factor
c=imcrop(a,[1 1 dim  dim]);
cr=c(:,:,1);
cg=c(:,:,2);
cb=c(:,:,3);
[count,x2r] = imhist(cr);
countcr = count;
[count,x2g] = imhist(cg);
countcg = count;
[count,x2b] = imhist(cb);
countcb = count;
countclr = (1/mx2r)*countcr;
countclg = (1/mx2g)*countcg;
countclb = (1/mx2b)*countcb;
convimr=zeros(diml,diml);
convimg=zeros(diml,diml);
convimb=zeros(diml,diml);

%membership value of red, green, and blue channels in the image windows%
for i1=1:dim:diml
    for i2=1:dim:diml
        dr=[];dg=[];db=[];
        countdivr=0.0;
        countdivg=0.0;
        countdivb=0.0;
        im = imcrop(a,[i2 i1 dim-1 dim-1]);
        imr=im(:,:,1);
        img=im(:,:,2);
        imb=im(:,:,3);
        imlr=double(imr);
        imlg=double(img);
        imlb=double(imb);
        [count,xlr]=imhist(imr);
        countimr=count;
        [count,xlg]=imhist(img);
        countimg=count;
        [count,xlb]=imhist(imb);
        countimb=count;
        countimlr=(1/mx2r)*countimr;
        countimlg=(1/mx2g)*countimg;
        countimlb=(1/mx2b)*countimb;
       % computing the fuzzy ratio between the each database image and
       % query image for each channel using fuzzy similarity
       for p=1:256
           if countimr(p)>0 ||countcr(p) >0
           counr=(min(countclr(p),countimlr(p)))./(min(countclr(p),countimlr(p))+0.5*min(countclr(p),1- ...
           countimlr(p))+0.5*min(countimlr(p),1-countclr(p)));
           countdivr=[dr;counr];
           dr=countdivr;
           else
               countdivr=0;
               countdivr=[dr;countdivr];
               dr=countdivr;
           end
           if countimg(p)>0 | countcg(p)>0 
           coung=(min(countclg(p),countimlg(p)))./(min(countclg(p),countimlg(p))+0.5*min(countclg(p),1- ...
           countimlg(p))+0.5*min(countimlg(p),1-countclg(p)));
           countdivg=[dg;coung];
           dg=countdivg;
           else
               countdivg=0;
               countdivg=[dg;countdivg];
               dg=countdivg;
           end
           if countimb(p)>0 | countcb(p)>0 
           counb=(min(countclb(p),countimlb(p)))./(min(countclb(p),countimlb(p))+0.5*min(countclb(p),1- ...
           countimlb(p))+0.5*min(countimlb(p),1-countclb(p)));
           countdivb=[db;counb];
           db=countdivb;
           else
               countdivb=0;
               countdivb=[db;countdivb];
               db=countdivb;
           end
       end
       hr = (countdivr);
       hg = (countdivg);
       hb = (countdivb);
       %back propogation of fuzzy ratio on the images%
       for i=i1:i1+(dim-1)
           for j=i2:i2+(dim-1)
               sr=imlr(i-i1+1,j-i2+1);
               convimr(i,j)=hr(sr+1);
               sg=imlg(i-i1+1,j-i2+1);
               convimg(i,j)=hg(sg+1);
               sb=imlb(i-i1+1,j-i2+1);
               convimb(i,j)=hb(sb+1);
           end
       end
    end
end
convimr;convimg;convimb;
tempr=max(max(convimr));
tempg=max(max(convimg));
tempb=max(max(convimb));
finalimr = (1/tempr)*255*convimr;finalimrl=uint8(finalimr);
finalimg = (1/tempg)*255*convimg;finalimgl=uint8(finalimg);
finalimb = (1/tempb)*255*convimb;finalimbl=uint8(finalimb);
[count,x9r]=imhist(finalimrl);

%finding the highest peak for each channel%
countxyr=count;
mxconr=max(countxyr);
[x9r]=find(mxconr==countxyr);
[count,x9g]=imhist(finalimgl);
countxyg=count;
mxcong=max(countxyg);
[x9g]=find(mxcong==countxyg);
[count,x9b]=imhist(finalimbl);
countxyb=count;
mxconb=max(countxyb);
[x9b]=find(mxconb==countxyb);

%thresholding at the highest peak for each channel%
for w1 = 1:diml
    for w2 = 1:diml
        if finalimrl(w1,w2) == (x9r-1)
            segimr(w1,w2) = 0;
        else segimr(w1,w2) = 1;
        end
    end
end
for w1 = 1:diml
    for w2 = 1:diml
        if finalimgl(w1,w2) == (x9r-1)
            segimg(w1,w2) = 0;
        else segimg(w1,w2) = 1;
        end
    end
end
for w1 = 1:diml
    for w2 = 1:diml
        if finalimbl(w1,w2) == (x9r-1)
            segimb(w1,w2) = 0;
        else segimb(w1,w2) = 1;
        end
    end
end

%figure, imshow(b)
figure, imshow(c)
image1=segimr + segimg+ segimb;
tmp = max(max(image1));
image = (1/tmp)*image1;
[count,x10] = imhist(image);
countim = count;
mxc = max(countim);
[x10] = find(mxc == countim);
for w3 = 1:diml
    for w4 = 1:diml
        if image(w3,w4) == (x10-1)/255
            finish(w3,w4) = 0;
        else finim(w3,w4) = 1;
        end
    end
end

% filtering with Gaussian filter%
h= fspecial('gaussian',3,4);
finalimage = filter2(h,finim);
figure, imshow(finalimage)

toc
