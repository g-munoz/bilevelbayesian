Fx = [ 16.0; -47.0; 13.0; 2.0; 39.0];
Fy = [ -28.0; -17.0; 3.0; -49.0; -15.0; 13.0; 5.0; -42.0; 48.0; 16.0; 13.0; -41.0; -31.0; -26.0; -31.0];
Gx = [];
Gy = [];
gx = [ 18.0 0.0 0.0 6.0 27.0; 12.0 19.0 8.0 0.0 35.0; 27.0 27.0 26.0 0.0 40.0; 0.0 5.0 49.0 37.0 8.0; 28.0 37.0 45.0 6.0 11.0; 32.0 46.0 10.0 49.0 39.0; 0.0 0.0 0.0 0.0 45.0; 0.0 18.0 0.0 39.0 10.0; 0.0 0.0 47.0 42.0 0.0; 0.0 28.0 45.0 36.0 0.0; 40.0 32.0 2.0 21.0 0.0; 0.0 0.0 47.0 0.0 15.0; 32.0 0.0 27.0 0.0 48.0; 0.0 50.0 33.0 28.0 22.0; 0.0 46.0 12.0 0.0 30.0; 6.0 0.0 50.0 13.0 41.0; 24.0 0.0 20.0 1.0 33.0; 0.0 0.0 40.0 0.0 13.0; 0.0 42.0 13.0 26.0 35.0; 0.0 22.0 28.0 38.0 6.0; 1.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 -1.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0];
gy = [ 2.0 49.0 11.0 34.0 0.0 0.0 0.0 5.0 26.0 0.0 33.0 50.0 47.0 24.0 0.0; 12.0 30.0 8.0 2.0 0.0 2.0 34.0 2.0 36.0 0.0 42.0 31.0 37.0 11.0 42.0; 10.0 0.0 0.0 26.0 0.0 0.0 49.0 0.0 0.0 45.0 24.0 41.0 0.0 39.0 31.0; 5.0 0.0 37.0 47.0 0.0 31.0 43.0 23.0 41.0 12.0 33.0 44.0 34.0 49.0 6.0; 0.0 0.0 3.0 8.0 0.0 9.0 29.0 16.0 24.0 50.0 43.0 3.0 44.0 6.0 25.0; 0.0 38.0 0.0 38.0 4.0 39.0 37.0 0.0 0.0 0.0 0.0 0.0 14.0 8.0 22.0; 0.0 42.0 0.0 34.0 37.0 0.0 0.0 0.0 0.0 31.0 0.0 0.0 20.0 35.0 29.0; 0.0 34.0 13.0 0.0 18.0 4.0 45.0 0.0 39.0 39.0 17.0 13.0 5.0 31.0 0.0; 0.0 39.0 0.0 18.0 0.0 28.0 36.0 36.0 41.0 0.0 0.0 16.0 0.0 1.0 36.0; 28.0 0.0 0.0 14.0 14.0 0.0 35.0 5.0 30.0 0.0 41.0 18.0 5.0 23.0 0.0; 22.0 0.0 45.0 0.0 0.0 18.0 35.0 0.0 0.0 0.0 23.0 8.0 0.0 48.0 11.0; 24.0 24.0 0.0 0.0 0.0 0.0 0.0 42.0 9.0 0.0 17.0 15.0 0.0 18.0 27.0; 13.0 12.0 19.0 0.0 31.0 3.0 49.0 1.0 8.0 25.0 37.0 9.0 18.0 18.0 0.0; 45.0 23.0 43.0 15.0 20.0 35.0 31.0 23.0 49.0 15.0 26.0 7.0 25.0 44.0 49.0; 2.0 0.0 0.0 47.0 17.0 0.0 0.0 0.0 0.0 32.0 32.0 0.0 32.0 10.0 23.0; 20.0 30.0 20.0 13.0 0.0 39.0 41.0 5.0 12.0 11.0 44.0 23.0 14.0 45.0 8.0; 41.0 38.0 30.0 0.0 0.0 8.0 43.0 31.0 22.0 1.0 18.0 42.0 0.0 0.0 42.0; 0.0 31.0 0.0 0.0 17.0 3.0 41.0 18.0 21.0 33.0 1.0 0.0 0.0 0.0 9.0; 36.0 1.0 44.0 46.0 0.0 0.0 25.0 0.0 27.0 0.0 43.0 44.0 14.0 42.0 14.0; 44.0 44.0 8.0 34.0 28.0 0.0 10.0 2.0 25.0 9.0 15.0 5.0 5.0 0.0 46.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0];
bg = [ -770.0; -667.0; -703.0; -741.0; -722.0; -531.0; -240.0; -359.0; -112.0; -470.0; -679.0; -198.0; -481.0; -699.0; -526.0; -566.0; -392.0; -6.0; -797.0; -205.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0];
bG = [];
