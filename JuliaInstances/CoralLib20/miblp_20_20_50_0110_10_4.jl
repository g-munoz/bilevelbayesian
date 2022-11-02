Fx = [ -40.0; -35.0; 43.0; 12.0; -12.0; -16.0; 10.0; -4.0; -27.0; 13.0];
Fy = [ -26.0; 20.0; 21.0; -49.0; -13.0; -2.0; -34.0; -30.0; -27.0; 23.0];
Gx = [];
Gy = [];
gx = [ 0.0 0.0 0.0 0.0 45.0 0.0 0.0 0.0 4.0 0.0; 11.0 0.0 21.0 25.0 39.0 48.0 27.0 32.0 28.0 20.0; 13.0 0.0 17.0 50.0 0.0 0.0 31.0 24.0 0.0 0.0; 10.0 10.0 4.0 16.0 16.0 0.0 0.0 34.0 7.0 45.0; 5.0 10.0 32.0 2.0 34.0 22.0 37.0 50.0 0.0 13.0; 0.0 5.0 0.0 0.0 47.0 0.0 48.0 19.0 19.0 0.0; 23.0 50.0 0.0 17.0 0.0 17.0 32.0 40.0 49.0 41.0; 1.0 14.0 40.0 38.0 3.0 8.0 24.0 17.0 26.0 0.0; 0.0 32.0 19.0 23.0 0.0 26.0 35.0 36.0 25.0 17.0; 0.0 7.0 0.0 32.0 50.0 29.0 34.0 0.0 0.0 13.0; 49.0 44.0 43.0 15.0 0.0 8.0 0.0 0.0 2.0 0.0; 3.0 35.0 30.0 28.0 0.0 44.0 22.0 49.0 11.0 24.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 32.0 26.0 16.0; 0.0 3.0 30.0 0.0 11.0 11.0 0.0 7.0 0.0 0.0; 44.0 0.0 0.0 19.0 47.0 0.0 0.0 21.0 36.0 41.0; 0.0 34.0 5.0 6.0 47.0 25.0 46.0 38.0 46.0 39.0; 5.0 25.0 40.0 45.0 0.0 47.0 47.0 27.0 0.0 3.0; 2.0 25.0 36.0 21.0 47.0 22.0 16.0 0.0 0.0 25.0; 0.0 32.0 0.0 33.0 48.0 0.0 0.0 0.0 2.0 26.0; 46.0 40.0 40.0 37.0 46.0 13.0 37.0 0.0 28.0 18.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
gy = [ 7.0 2.0 0.0 0.0 0.0 38.0 27.0 21.0 0.0 44.0; 36.0 27.0 37.0 42.0 16.0 41.0 34.0 34.0 42.0 12.0; 0.0 41.0 0.0 0.0 32.0 15.0 0.0 0.0 0.0 35.0; 44.0 0.0 27.0 0.0 0.0 37.0 33.0 0.0 10.0 0.0; 12.0 8.0 30.0 43.0 42.0 2.0 1.0 23.0 49.0 27.0; 27.0 17.0 38.0 0.0 0.0 0.0 35.0 33.0 0.0 9.0; 44.0 0.0 0.0 35.0 0.0 13.0 40.0 30.0 10.0 34.0; 45.0 17.0 20.0 0.0 0.0 2.0 9.0 0.0 0.0 12.0; 45.0 28.0 42.0 48.0 13.0 24.0 24.0 30.0 14.0 0.0; 4.0 0.0 0.0 0.0 32.0 8.0 0.0 9.0 50.0 0.0; 35.0 17.0 29.0 16.0 0.0 14.0 37.0 0.0 27.0 0.0; 1.0 47.0 34.0 8.0 20.0 10.0 5.0 8.0 0.0 47.0; 48.0 0.0 35.0 38.0 27.0 4.0 9.0 0.0 0.0 10.0; 0.0 0.0 0.0 6.0 31.0 0.0 21.0 0.0 0.0 28.0; 47.0 14.0 46.0 24.0 0.0 0.0 30.0 22.0 9.0 6.0; 30.0 49.0 27.0 44.0 16.0 0.0 0.0 44.0 37.0 38.0; 44.0 33.0 13.0 50.0 36.0 10.0 30.0 35.0 44.0 0.0; 0.0 13.0 39.0 0.0 30.0 40.0 11.0 23.0 29.0 33.0; 2.0 11.0 41.0 0.0 0.0 0.0 5.0 0.0 0.0 27.0; 0.0 6.0 37.0 0.0 40.0 16.0 0.0 8.0 26.0 0.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0];
bg = [ -218.0; -675.0; -268.0; -410.0; -802.0; -317.0; -826.0; -218.0; -807.0; -246.0; -445.0; -789.0; -593.0; -281.0; -623.0; -1024.0; -669.0; -743.0; -603.0; -651.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0];
bG = [];