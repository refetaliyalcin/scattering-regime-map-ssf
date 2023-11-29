function [S] = SSF_correction(f_v, theta, lamda, r)

% Source: https://bitbucket.org/planetarysystemresearch/rtcb_public/src/master/src/mie.f90
% #####################################################################
% Copyright (c) 2015,2016,2017 Karri Muinonen, Timo Väisänen, 
% Antti Penttilä and University of Helsinki
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above
%       copyright notice, this list of conditions and the following
%       disclaimer in the documentation and/or other materials provided
%       with the distribution.
%     * Neither the name of University of Helsinki nor the names of
%       its contributors may be used to endorse or promote products
%       derived from this software without specific prior written
%       permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


    S=zeros(length(theta),1);
    for i1 = 1: size(theta)
        q = 4 * pi * sin(theta(i1)/2) / lamda;
        x = q * r;
        psi = sin(x) / x;
        sigma = 3*(sin(x)/x^3 - cos(x)/x^2);
        alpha = f_v/(1-f_v) * ((1 + 3*f_v/(1-f_v)) * sigma + 3*psi) + cos(x);
        beta = f_v/(1-f_v) * x * sigma + sin(x);
        S(i1)=1 / (alpha^2 + beta^2);

        % %different formula
        % % n = f_v *3/4 / pi /r^3;
        % 
        % p = 4*pi*sin(theta(i1)/2) / lamda;
        % 
        % alpha = (1 + 2*f_v)^2 / (1 - f_v)^4;
        % beta = -6*f_v*(1 + f_v/2)^2 / (1 - f_v)^4;
        % delta = alpha*f_v/2;
        % %delta = f_v*(1 + 2*f_v)^2 /(2*(1 - f_v)^2)
        % u=2*p*r;
        % 
        % if(abs(theta(i1)) <  10^-8)
        %     C = -24 * f_v * (alpha/3 + beta/4 + delta/6);
        % else
        %     C = 24 *(f_v) * ((alpha+beta+delta) * cos(u)/u^2 -(alpha+ 2*beta + 4*delta) * sin(u)/u^3 -  2*(beta+6*delta)*cos(u)/u^4 + 2*beta/u^4 + 24*delta*sin(u)/u^5 + 24*delta*(cos(u)-1)/u^6);
        % end
        % 
        % S(i1) = 1/(1-C);

    end
end