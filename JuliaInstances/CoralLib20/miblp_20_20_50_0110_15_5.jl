Fx = [ 4.0; -3.0; 18.0; 44.0; 39.0];
Fy = [ 16.0; 9.0; 41.0; -10.0; 42.0; -26.0; -41.0; 43.0; 32.0; 27.0; 42.0; -37.0; -21.0; 37.0; 6.0];
Gx = [];
Gy = [];
gx = [ 0.0 11.0 23.0 47.0 7.0; 32.0 0.0 3.0 0.0 13.0; 0.0 0.0 50.0 0.0 7.0; 37.0 0.0 21.0 19.0 35.0; 34.0 14.0 8.0 37.0 20.0; 14.0 33.0 46.0 44.0 12.0; 36.0 36.0 2.0 49.0 18.0; 3.0 24.0 28.0 1.0 39.0; 0.0 0.0 0.0 46.0 2.0; 0.0 7.0 9.0 0.0 0.0; 0.0 0.0 0.0 19.0 38.0; 23.0 0.0 1.0 9.0 1.0; 0.0 34.0 38.0 43.0 0.0; 22.0 0.0 49.0 0.0 19.0; 0.0 3.0 34.0 49.0 30.0; 40.0 28.0 40.0 0.0 37.0; 6.0 45.0 50.0 10.0 34.0; 0.0 3.0 43.0 45.0 5.0; 0.0 18.0 0.0 31.0 40.0; 35.0 21.0 0.0 31.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 -1.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0];
gy = [ 26.0 10.0 44.0 34.0 43.0 38.0 45.0 32.0 41.0 48.0 11.0 0.0 43.0 2.0 0.0; 41.0 16.0 48.0 0.0 0.0 0.0 0.0 16.0 8.0 41.0 0.0 0.0 0.0 19.0 0.0; 0.0 13.0 0.0 30.0 0.0 29.0 11.0 46.0 0.0 38.0 0.0 1.0 0.0 0.0 0.0; 14.0 18.0 11.0 3.0 0.0 48.0 0.0 35.0 3.0 2.0 8.0 29.0 26.0 42.0 10.0; 50.0 25.0 44.0 7.0 15.0 29.0 29.0 0.0 9.0 21.0 8.0 29.0 8.0 0.0 11.0; 37.0 0.0 0.0 47.0 32.0 0.0 0.0 4.0 11.0 0.0 8.0 49.0 30.0 23.0 48.0; 48.0 49.0 0.0 47.0 47.0 0.0 49.0 41.0 36.0 33.0 24.0 0.0 0.0 46.0 31.0; 10.0 38.0 40.0 49.0 47.0 20.0 5.0 6.0 47.0 3.0 38.0 46.0 36.0 50.0 7.0; 23.0 41.0 46.0 5.0 1.0 0.0 0.0 8.0 33.0 43.0 0.0 10.0 30.0 13.0 14.0; 23.0 39.0 0.0 0.0 48.0 10.0 8.0 49.0 11.0 8.0 6.0 29.0 0.0 32.0 0.0; 0.0 0.0 0.0 0.0 0.0 3.0 32.0 0.0 10.0 0.0 35.0 30.0 0.0 33.0 0.0; 0.0 0.0 19.0 14.0 0.0 0.0 46.0 0.0 7.0 0.0 9.0 30.0 14.0 26.0 35.0; 0.0 19.0 0.0 36.0 8.0 30.0 28.0 11.0 26.0 0.0 7.0 0.0 10.0 26.0 0.0; 0.0 0.0 0.0 1.0 0.0 26.0 22.0 0.0 31.0 0.0 0.0 3.0 0.0 0.0 0.0; 0.0 47.0 0.0 0.0 0.0 0.0 16.0 49.0 30.0 39.0 30.0 0.0 38.0 0.0 12.0; 18.0 9.0 0.0 13.0 37.0 13.0 33.0 4.0 50.0 0.0 18.0 15.0 5.0 38.0 26.0; 4.0 8.0 33.0 0.0 11.0 34.0 34.0 2.0 2.0 37.0 33.0 0.0 39.0 8.0 49.0; 16.0 8.0 0.0 49.0 0.0 37.0 15.0 0.0 36.0 17.0 19.0 9.0 0.0 7.0 28.0; 23.0 0.0 0.0 34.0 0.0 28.0 0.0 0.0 47.0 50.0 36.0 45.0 8.0 15.0 35.0; 0.0 9.0 14.0 37.0 31.0 36.0 27.0 28.0 7.0 5.0 39.0 0.0 33.0 25.0 38.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0];
bg = [ -931.0; -216.0; -612.0; -649.0; -763.0; -955.0; -1162.0; -1038.0; -481.0; -673.0; -535.0; -599.0; -653.0; -413.0; -876.0; -839.0; -964.0; -799.0; -839.0; -807.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0];
bG = [];