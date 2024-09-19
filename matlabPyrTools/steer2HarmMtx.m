% MTX = steer2HarmMtx(HARMONICS, ANGLES, REL_PHASES)
%
% Compute a steering matrix (maps a directional basis set onto the
% angular Fourier harmonics).  
% 计算一个转向矩阵(将一个方向基映射到角傅里叶谐波上)。
% HARMONICS is a vector specifying the
% angular harmonics contained in the steerable basis/filters.  
% HARMONICS是一个矢量，指定可操纵基/滤波器中包含的角谐波。
% ANGLES 
% (optional) is a vector specifying the angular position of each filter. 
% ANGLES(可选)是指定每个过滤器的角度位置的矢量。
% REL_PHASES (optional, default = 'even') specifies whether the harmonics
% REL_PHASES(可选，default = 'even')指定谐波是沿这些位置对齐的余弦相位还是正弦相位。
% are cosine or sine phase aligned about those positions.
% The result matrix is suitable for passing to the function STEER.

% Eero Simoncelli, 7/96.

function mtx = steer2HarmMtx(harmonics, angles, evenorodd)

%%=================================================================
%%% Optional Parameters:

if (exist('evenorodd') ~= 1)
  evenorodd = 'even';
end

% Make HARMONICS a row vector
harmonics = harmonics(:)';

numh = 2*size(harmonics,2) - any(harmonics == 0);

if (exist('angles') ~= 1)
  angles = pi * [0:numh-1]'/numh;
else
  angles = angles(:);
end

%%=================================================================

if isstr(evenorodd)
  if strcmp(evenorodd,'even')
    evenorodd = 0;
  elseif strcmp(evenorodd,'odd')
    evenorodd = 1;
  else
    error('EVEN_OR_ODD should be the string  EVEN or ODD');
  end
end

%% Compute inverse matrix, which maps Fourier components onto 
%% steerable basis.
imtx = zeros(size(angles,1),numh);
col = 1;
for h=harmonics
  args = h*angles;
  if (h == 0)
    imtx(:,col) = ones(size(angles));
    col = col+1;
  elseif evenorodd
    imtx(:,col) = sin(args);
    imtx(:,col+1) = -cos(args);
    col = col+2;
  else
    imtx(:,col) = cos(args);
    imtx(:,col+1) = sin(args);
    col = col+2;
  end
end
  
r = rank(imtx);
if (( r ~= numh ) & ( r ~= size(angles,1) ))
  fprintf(2,'WARNING: matrix is not full rank');
end  

mtx = pinv(imtx);

