% Series of examples to demonstrate the use of ParticleSpherical.
%
% See also Particle, ParticleSpherical.
%
% The OTGO - Optical Tweezers in Geometrical Optics
% software package complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

%   Author: Giovanni Volpe
%   Version: 1.0.0
%   Date: 2014/01/01


example('Use of ParticleSpherical')

%% DEFINITION OF PARTICLESPHERICAL
exampletitle('DEFINITION OF PARTICLESPHERICAL')

examplecode('c = Point(1,1,1);')
examplecode('r = 1;')
examplecode('nm = 1;')
examplecode('np = 1.5;')
examplecode('bead = ParticleSpherical(c,r,nm,np)')
examplewait()

%% PLOTTING OF PARTICLESPHERICAL
exampletitle('PLOTTING OF PARTICLESPHERICAL')

figure
title('PARTICLESPHERICAL')
hold on
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')

examplecode('bead.plot();')
examplewait()

%% SCATTERING
exampletitle('SCATTERING')

examplecode('mr = 3;')
examplecode('nr = 2;')
examplecode('v = Vector(zeros(mr,nr),zeros(mr,nr),zeros(mr,nr),rand(mr,nr),rand(mr,nr),rand(mr,nr));')
examplecode('P = ones(mr,nr);')
examplecode('pol = Vector(zeros(mr,nr),zeros(mr,nr),zeros(mr,nr),ones(mr,nr),ones(mr,nr),ones(mr,nr)); pol = v*pol;')
examplecode('r = Ray(v,P,pol);')
examplewait()

examplecode('r.plot(''color'',''k'');')
examplewait()

examplecode('r_vec = bead.scattering(r)')
examplewait()

examplecode('rr = r_vec(1).r;')
examplecode('rr.plot(''color'',''r'');')
examplecode('rt = r_vec(1).t;')
examplecode('rt.plot(''color'',''b'');')
examplewait()

examplecode('rr = r_vec(2).r;')
examplecode('rr.plot(''color'',''r'');')
examplecode('rt = r_vec(2).t;')
examplecode('rt.plot(''color'',''b'');')
examplewait()

examplecode('rr = r_vec(3).r;')
examplecode('rr.plot(''color'',''r'');')
examplecode('rt = r_vec(3).t;')
examplecode('rt.plot(''color'',''b'');')
examplewait()

examplecode('rr = r_vec(4).r;')
examplecode('rr.plot(''color'',''r'');')
examplecode('rt = r_vec(4).t;')
examplecode('rt.plot(''color'',''b'');')
examplewait()

examplecode('rr = r_vec(5).r;')
examplecode('rr.plot(''color'',''r'');')
examplecode('rt = r_vec(5).t;')
examplecode('rt.plot(''color'',''b'');')
examplewait()

%% FORCE
exampletitle('FORCE')

examplecode('F = bead.force(r)*1e+15 % fN')
examplewait()

%% TORQUE
exampletitle('TORQUE')

examplecode('T = bead.torque(r,1e-21)*1e+21 % fN*nm')