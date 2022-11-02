Fx = [ 3.0; -47.0; -38.0; -13.0; -39.0; 35.0; 48.0; 32.0; 46.0; 6.0];
Fy = [ 7.0; 30.0; -6.0; 8.0; -24.0; 9.0; 44.0; 24.0; -36.0; 19.0];
Gx = [];
Gy = [];
gx = [ 24.0 20.0 0.0 0.0 8.0 47.0 43.0 8.0 0.0 0.0; 13.0 33.0 0.0 0.0 0.0 2.0 33.0 24.0 19.0 9.0; 0.0 0.0 0.0 0.0 0.0 23.0 4.0 0.0 0.0 0.0; 2.0 5.0 20.0 27.0 0.0 26.0 26.0 13.0 23.0 15.0; 0.0 0.0 0.0 27.0 47.0 2.0 13.0 16.0 36.0 0.0; 0.0 6.0 0.0 7.0 0.0 20.0 23.0 36.0 0.0 0.0; 26.0 39.0 19.0 39.0 12.0 22.0 39.0 16.0 28.0 18.0; 6.0 0.0 2.0 17.0 32.0 27.0 44.0 0.0 0.0 20.0; 14.0 19.0 27.0 12.0 9.0 0.0 3.0 10.0 21.0 0.0; 25.0 0.0 0.0 0.0 2.0 23.0 0.0 34.0 0.0 0.0; 37.0 15.0 44.0 50.0 16.0 44.0 13.0 45.0 31.0 2.0; 7.0 0.0 0.0 33.0 0.0 1.0 0.0 50.0 0.0 24.0; 27.0 5.0 43.0 5.0 0.0 15.0 0.0 0.0 29.0 7.0; 17.0 0.0 46.0 0.0 43.0 18.0 40.0 3.0 0.0 17.0; 36.0 25.0 23.0 0.0 50.0 0.0 30.0 16.0 0.0 0.0; 37.0 6.0 0.0 33.0 20.0 0.0 0.0 36.0 0.0 30.0; 0.0 5.0 13.0 27.0 0.0 49.0 0.0 18.0 0.0 8.0; 17.0 17.0 7.0 41.0 16.0 0.0 0.0 37.0 35.0 3.0; 17.0 0.0 10.0 21.0 0.0 36.0 47.0 0.0 0.0 25.0; 23.0 38.0 33.0 9.0 12.0 45.0 17.0 45.0 42.0 19.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
gy = [ 28.0 13.0 0.0 36.0 0.0 33.0 0.0 0.0 0.0 0.0; 19.0 0.0 0.0 0.0 24.0 49.0 39.0 0.0 0.0 46.0; 40.0 3.0 21.0 2.0 0.0 0.0 5.0 43.0 0.0 43.0; 16.0 50.0 45.0 0.0 26.0 30.0 0.0 33.0 17.0 16.0; 0.0 41.0 0.0 19.0 39.0 15.0 0.0 0.0 45.0 30.0; 18.0 0.0 0.0 0.0 21.0 6.0 0.0 18.0 16.0 0.0; 3.0 24.0 15.0 9.0 9.0 44.0 4.0 50.0 42.0 10.0; 28.0 2.0 30.0 30.0 14.0 0.0 46.0 18.0 31.0 22.0; 7.0 50.0 49.0 0.0 32.0 2.0 47.0 27.0 24.0 49.0; 26.0 24.0 0.0 25.0 25.0 0.0 26.0 0.0 15.0 0.0; 3.0 9.0 46.0 20.0 24.0 1.0 22.0 6.0 20.0 41.0; 0.0 0.0 28.0 0.0 0.0 40.0 22.0 0.0 0.0 0.0; 0.0 5.0 41.0 39.0 28.0 0.0 0.0 0.0 19.0 3.0; 30.0 49.0 35.0 39.0 47.0 21.0 11.0 38.0 30.0 18.0; 3.0 0.0 21.0 23.0 0.0 41.0 46.0 50.0 0.0 49.0; 0.0 20.0 0.0 0.0 4.0 0.0 0.0 0.0 41.0 26.0; 30.0 0.0 7.0 0.0 9.0 39.0 25.0 0.0 1.0 25.0; 15.0 40.0 46.0 37.0 8.0 0.0 0.0 16.0 4.0 2.0; 49.0 19.0 0.0 3.0 5.0 28.0 19.0 13.0 50.0 0.0; 34.0 18.0 36.0 7.0 43.0 32.0 12.0 10.0 35.0 29.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0];
bg = [ -751.0; -427.0; -514.0; -747.0; -440.0; -436.0; -1061.0; -822.0; -399.0; -275.0; -832.0; -234.0; -315.0; -937.0; -770.0; -273.0; -588.0; -570.0; -860.0; -864.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0];
bG = [];
