function V = NumSoln(nx, ny, c, BcLeft, BcRight, BcTop, BcBottom)
%using the finite difference method


%Tariq Aboushaer
%101064544

    global C;
    g = sparse(nx*ny, ny*nx);
    b = zeros(1, nx*ny);
    for x=1:nx
        for y=1:ny
            n = y + (x - 1)*ny;
            nxm = y + (x - 2)*ny;
            nxp = y + (x)*ny;
            nym = (y-1) + (x - 1)*ny;
            nyp = (y+1) + (x - 1)*ny;

            if (x == 1 && y == 1)
                if (BcLeft == Inf)
                    rxp = (c(x,y) + c(x+1,y))/2;
                    ryp = (c(x,y) + c(x,y+1))/2;

                    g(n,n)   = -(rxp+ryp);
                    g(n,nxp) =  rxp;
                    g(n,nyp) =  ryp;
                else
                    g(n,n) = 1;
                    b(n) = BcLeft;
                end
            elseif (x == 1 && y == ny)
                if (BcLeft == Inf)
                    rxp = (c(x,y) + c(x+1,y))/2;
                    rym = (c(x,y) + c(x,y-1))/2;

                    g(n,n)   = -(rxp+rym);
                    g(n,nxp) =  rxp;
                    g(n,nym) =  rym;
                else
                    g(n,n) = 1;
                    b(n) = BcLeft;
                end
            elseif x == nx && y == 1 % Right side
                if (BcRight == Inf)
                    rxm = (c(x,y) + c(x-1,y))/2;
                    ryp = (c(x,y) + c(x,y+1))/2;
                    g(n,n)   = -(rxm+ryp);
                    g(n,nxm) =  rxm;
                    g(n,nyp) =  ryp;
                else
                    g(n,n) = 1;
                    b(n) = BcRight;
                end
            elseif x == nx && y == ny % Right side
                if (BcRight == Inf)
                    rxm = (c(x,y) + c(x-1,y))/2;
                    rym = (c(x,y) + c(x,y-1))/2;
                    g(n,n)   = -(rxm+rym);
                    g(n,nxm) =  rxm;
                    g(n,nym) =  rym;
                else
                    g(n,n) = 1;
                    b(n) = BcRight;
                end
            elseif (x == 1) % Left Side
                if (BcLeft == Inf)
                    rxp = (c(x,y) + c(x+1,y))/2;
                    rym = (c(x,y) + c(x,y-1))/2;
                    ryp = (c(x,y) + c(x,y+1))/2;

                    g(n,n)   = -(rxp+rym+ryp);
                    g(n,nxp) =  rxp;
                    g(n,nym) =  rym;
                    g(n,nyp) =  ryp;
                else
                    g(n,n) = 1;
                    b(n) = BcLeft;
                end
            elseif x == nx % Right side
                if (BcRight == Inf)
                    rxm = (c(x,y) + c(x-1,y))/2;
                    rym = (c(x,y) + c(x,y-1))/2;
                    ryp = (c(x,y) + c(x,y+1))/2;
                    g(n,n)   = -(rxm+rym+ryp);
                    g(n,nxm) =  rxm;
                    g(n,nym) =  rym;
                    g(n,nyp) =  ryp;
                else
                    g(n,n) = 1;
                    b(n) = BcRight;
                end
            elseif y == 1 % Top Side
                if (BcTop == Inf)
                    rxm = (c(x,y) + c(x-1,y))/2;
                    rxp = (c(x,y) + c(x+1,y))/2;
                    ryp = (c(x,y) + c(x,y+1))/2;
                    g(n,n) = -(rxm+rxp+ryp);
                    g(n,nxm) =  rxm;
                    g(n,nxp) =  rxp;
                    g(n,nyp) =  ryp;
                else
                    g(n,n) = 1;
                    b(n) = BcTop;
                end
            elseif y == ny % Bottom Side
                if (BcBottom == Inf)
                    rxm = (c(x,y) + c(x-1,y))/2;
                    rxp = (c(x,y) + c(x+1,y))/2;
                    rym = (c(x,y) + c(x,y-1))/2;
                    g(n,n) = -(rxm+rxp+rym);
                    g(n,nxm) =  rxm;
                    g(n,nxp) =  rxp;
                    g(n,nym) =  rym;
                else
                    g(n,n) = 1;
                    b(n) = BcBottom;
                end
            else % Bulk Area
                rxm = (c(x,y) + c(x-1,y))/2;
                rxp = (c(x,y) + c(x+1,y))/2;
                rym = (c(x,y) + c(x,y-1))/2;
                ryp = (c(x,y) + c(x,y+1))/2;
                
                g(n,n)   = -(rxm+rxp+rym+ryp);
                g(n,nxm) =  rxm;
                g(n,nxp) =  rxp;
                g(n,nym) =  rym;
                g(n,nyp) =  ryp;
            end
        end
    end
    
    V_temp = g\b';
    
    V = zeros(nx,ny,1);
    for x=1:nx
        for y=1:ny
            V(x,y) = V_temp(y + (x - 1)*ny);
        end
    end
end
