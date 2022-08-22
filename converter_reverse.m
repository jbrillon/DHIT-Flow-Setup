function converter_2(NP)
%MUCH FASTER THAN CONVERTER 1
for NUM = 8000:100:8000
    A = [];
    for i=0:NP-1
        horzcat(horzcat(horzcat(num2str(1000000+NUM),'__'),num2str(i)),'.dat_equidist.dat')
        B = importdata(horzcat(horzcat(horzcat(num2str(1000000+NUM),'__'),num2str(i)),'.dat_equidist.dat'));
        A = [A;B];
    end
    
    XU = unique(A(:,1));
    YU = unique(A(:,2));
    ZU = unique(A(:,3));

    disp('Number of Points')
    N  = max(size(XU));

    B  = zeros(N*N*N,min(size(A))+1);

    %Sort First Column
    %Pre-Matches the x Coordinates
    [Y,I] = sort(A(:,1));
    A     = A(I,:);
    
    %Count the number of values in A with different x's
    %1 -> x-coordinate
    %2 -> Number with it
    %3 -> Starting Location in A
    %4 -> Ending Location in A
    for i=1:max(size(XU))
        XU(i,2) = 0;
        XU(i,3) = 0;
        XU(i,4) = 0;
        for j=1:max(size(A))
            if(A(j,1)==XU(i,1))
                XU(i,2) = XU(i,2) + 1;
            end
        end
        
        if(i==1)
            XU(i,3) = 1;
        end
        if(i>1)
            XU(i,3) = XU(i-1,4)+1;
        end
        XU(i,4) = XU(i,3) + XU(i,2)-1;
    end

    XU
    for i=1:N
        i
    for j=1:N
    for k=1:N
        x = XU(i,1);
        y = YU(j);
        z = ZU(k);

        for m=XU(i,3):XU(i,4)
            if(A(m,1) == x)
            if(A(m,2) == y)
            if(A(m,3) == z)
                B((k-1)*N*N+(j-1)*N+i,1) = B((k-1)*N*N+(j-1)*N+i,1) + A(m,1);
                B((k-1)*N*N+(j-1)*N+i,2) = B((k-1)*N*N+(j-1)*N+i,2) + A(m,2);
                B((k-1)*N*N+(j-1)*N+i,3) = B((k-1)*N*N+(j-1)*N+i,3) + A(m,3);
                B((k-1)*N*N+(j-1)*N+i,4) = B((k-1)*N*N+(j-1)*N+i,4) + A(m,4);
                B((k-1)*N*N+(j-1)*N+i,5) = B((k-1)*N*N+(j-1)*N+i,5) + A(m,5);
                B((k-1)*N*N+(j-1)*N+i,6) = B((k-1)*N*N+(j-1)*N+i,6) + A(m,6);
                B((k-1)*N*N+(j-1)*N+i,7) = B((k-1)*N*N+(j-1)*N+i,7) + A(m,7);
                B((k-1)*N*N+(j-1)*N+i,8) = B((k-1)*N*N+(j-1)*N+i,8) + 1;
            end
            end
            end
        end

    end
    end
    end
    
    max(size(A))

    for(i=1:max(size(B)))
        B(i,1) = B(i,1)/B(i,8);
        B(i,2) = B(i,2)/B(i,8);
        B(i,3) = B(i,3)/B(i,8);
        B(i,4) = B(i,4)/B(i,8);
        B(i,5) = B(i,5)/B(i,8);
        B(i,6) = B(i,6)/B(i,8);
        B(i,7) = B(i,7)/B(i,8);
    end

    fid = fopen(horzcat('velocity',horzcat(num2str(1000000+NUM)),'.fld'),'w');
    for i=1:max(size(B))
        fprintf(fid,'%19.16f\n%19.16f\n%19.16f\n',B(i,5)/B(i,4),B(i,6)/B(i,4),B(i,7)/B(i,4));
    end
    fclose(fid);
    
end