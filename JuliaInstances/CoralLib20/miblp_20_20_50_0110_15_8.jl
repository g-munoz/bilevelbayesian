Fx = [ -18.0; 15.0; -34.0; -5.0; -9.0];
Fy = [ -46.0; -20.0; 11.0; -34.0; 49.0; 7.0; -6.0; -35.0; 31.0; -2.0; -44.0; 23.0; -19.0; -24.0; -10.0];
Gx = [];
Gy = [];
gx = [ 6.0 38.0 37.0 38.0 0.0; 47.0 27.0 44.0 26.0 27.0; 7.0 20.0 37.0 0.0 0.0; 0.0 41.0 5.0 4.0 35.0; 30.0 24.0 0.0 0.0 17.0; 14.0 37.0 38.0 9.0 45.0; 33.0 43.0 31.0 30.0 9.0; 0.0 0.0 0.0 14.0 0.0; 30.0 20.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 50.0 8.0 23.0 49.0 32.0; 0.0 0.0 23.0 10.0 33.0; 0.0 47.0 32.0 0.0 18.0; 8.0 14.0 0.0 50.0 39.0; 0.0 0.0 0.0 0.0 0.0; 41.0 0.0 33.0 24.0 36.0; 0.0 0.0 0.0 0.0 0.0; 37.0 36.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 48.0 11.0 32.0 0.0 42.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 -1.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0];
gy = [ 6.0 33.0 0.0 16.0 1.0 8.0 26.0 7.0 8.0 0.0 38.0 32.0 6.0 48.0 0.0; 22.0 21.0 37.0 0.0 17.0 11.0 0.0 38.0 12.0 10.0 15.0 0.0 21.0 14.0 43.0; 6.0 35.0 7.0 47.0 0.0 50.0 34.0 39.0 38.0 19.0 0.0 34.0 12.0 41.0 0.0; 20.0 0.0 46.0 37.0 11.0 0.0 0.0 0.0 42.0 21.0 0.0 34.0 19.0 13.0 24.0; 4.0 0.0 21.0 0.0 17.0 41.0 49.0 0.0 8.0 0.0 12.0 0.0 0.0 1.0 0.0; 37.0 19.0 33.0 10.0 43.0 47.0 9.0 30.0 20.0 22.0 15.0 0.0 28.0 16.0 32.0; 17.0 0.0 25.0 13.0 0.0 3.0 9.0 15.0 28.0 13.0 6.0 21.0 43.0 27.0 4.0; 0.0 6.0 0.0 44.0 0.0 0.0 37.0 0.0 23.0 0.0 0.0 0.0 23.0 0.0 0.0; 30.0 14.0 42.0 0.0 0.0 35.0 0.0 0.0 45.0 6.0 7.0 17.0 0.0 6.0 40.0; 17.0 0.0 35.0 11.0 0.0 2.0 29.0 0.0 0.0 0.0 18.0 10.0 47.0 0.0 25.0; 26.0 45.0 11.0 31.0 29.0 41.0 5.0 47.0 36.0 20.0 36.0 43.0 22.0 8.0 23.0; 0.0 0.0 24.0 16.0 24.0 0.0 0.0 0.0 0.0 44.0 0.0 48.0 0.0 39.0 14.0; 37.0 8.0 36.0 0.0 30.0 29.0 0.0 29.0 0.0 0.0 0.0 43.0 0.0 21.0 0.0; 5.0 50.0 45.0 0.0 0.0 33.0 8.0 37.0 21.0 35.0 6.0 7.0 21.0 41.0 19.0; 40.0 34.0 40.0 47.0 23.0 0.0 47.0 17.0 28.0 0.0 46.0 0.0 17.0 0.0 0.0; 35.0 45.0 0.0 30.0 11.0 18.0 37.0 24.0 0.0 0.0 0.0 0.0 44.0 0.0 34.0; 15.0 40.0 45.0 16.0 26.0 39.0 0.0 47.0 44.0 0.0 0.0 33.0 0.0 0.0 0.0; 43.0 20.0 3.0 30.0 26.0 32.0 23.0 0.0 9.0 12.0 7.0 0.0 0.0 0.0 0.0; 44.0 15.0 0.0 2.0 0.0 23.0 42.0 0.0 0.0 35.0 47.0 47.0 3.0 36.0 0.0; 13.0 0.0 0.0 4.0 0.0 7.0 5.0 29.0 11.0 0.0 0.0 14.0 2.0 1.0 22.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0];
bg = [ -810.0; -693.0; -790.0; -576.0; -457.0; -939.0; -764.0; -291.0; -674.0; -360.0; -1065.0; -238.0; -566.0; -806.0; -715.0; -596.0; -498.0; -589.0; -780.0; -281.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0];
bG = [];