Fx = [ 6.0; -29.0; 47.0; 24.0; 43.0; -37.0; -35.0; -15.0; -12.0; 28.0; 33.0; 20.0; -4.0; 40.0; -28.0];
Fy = [ 40.0; -36.0; 38.0; 31.0; -34.0];
Gx = [];
Gy = [];
gx = [ 26.0 41.0 42.0 1.0 25.0 6.0 24.0 0.0 43.0 28.0 31.0 39.0 41.0 50.0 30.0; 50.0 0.0 49.0 0.0 15.0 9.0 41.0 33.0 6.0 36.0 0.0 20.0 3.0 22.0 44.0; 20.0 35.0 50.0 0.0 0.0 0.0 33.0 36.0 30.0 47.0 7.0 0.0 35.0 40.0 9.0; 0.0 41.0 0.0 50.0 0.0 0.0 3.0 0.0 45.0 21.0 0.0 0.0 0.0 41.0 22.0; 32.0 1.0 36.0 0.0 46.0 0.0 5.0 0.0 19.0 0.0 0.0 0.0 36.0 12.0 42.0; 38.0 0.0 0.0 23.0 46.0 0.0 0.0 0.0 0.0 35.0 0.0 43.0 0.0 0.0 47.0; 0.0 38.0 14.0 17.0 9.0 0.0 10.0 1.0 0.0 1.0 0.0 48.0 0.0 36.0 0.0; 27.0 33.0 37.0 34.0 0.0 0.0 0.0 34.0 29.0 32.0 19.0 25.0 19.0 0.0 0.0; 0.0 50.0 20.0 34.0 36.0 0.0 2.0 0.0 0.0 2.0 23.0 21.0 0.0 0.0 7.0; 0.0 0.0 9.0 0.0 0.0 0.0 0.0 0.0 24.0 0.0 0.0 0.0 23.0 0.0 1.0; 1.0 15.0 24.0 49.0 0.0 35.0 41.0 0.0 0.0 0.0 0.0 0.0 18.0 22.0 0.0; 0.0 1.0 0.0 18.0 45.0 37.0 2.0 10.0 22.0 0.0 44.0 49.0 39.0 39.0 0.0; 15.0 44.0 48.0 35.0 44.0 0.0 13.0 34.0 0.0 0.0 30.0 0.0 10.0 13.0 45.0; 26.0 0.0 23.0 48.0 4.0 0.0 49.0 44.0 27.0 29.0 19.0 30.0 14.0 35.0 40.0; 27.0 40.0 49.0 41.0 20.0 7.0 15.0 0.0 0.0 0.0 28.0 31.0 0.0 21.0 41.0; 33.0 45.0 25.0 50.0 15.0 48.0 17.0 23.0 47.0 6.0 7.0 42.0 49.0 26.0 12.0; 2.0 33.0 0.0 0.0 0.0 0.0 25.0 30.0 0.0 16.0 48.0 0.0 5.0 0.0 29.0; 47.0 18.0 47.0 44.0 47.0 0.0 13.0 0.0 41.0 1.0 0.0 33.0 0.0 37.0 3.0; 0.0 0.0 0.0 0.0 32.0 13.0 28.0 24.0 46.0 15.0 44.0 0.0 0.0 50.0 0.0; 0.0 26.0 0.0 4.0 0.0 31.0 0.0 0.0 0.0 20.0 25.0 0.0 0.0 50.0 23.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0];
gy = [ 27.0 1.0 28.0 42.0 7.0; 3.0 0.0 19.0 2.0 43.0; 38.0 0.0 22.0 0.0 3.0; 0.0 33.0 45.0 34.0 6.0; 34.0 0.0 0.0 0.0 0.0; 0.0 0.0 29.0 0.0 32.0; 0.0 36.0 0.0 42.0 8.0; 32.0 0.0 0.0 37.0 0.0; 0.0 7.0 0.0 27.0 0.0; 14.0 22.0 0.0 11.0 0.0; 17.0 24.0 17.0 50.0 45.0; 5.0 45.0 48.0 6.0 16.0; 6.0 27.0 44.0 27.0 0.0; 21.0 38.0 39.0 26.0 0.0; 11.0 26.0 29.0 0.0 0.0; 40.0 39.0 8.0 49.0 42.0; 31.0 0.0 0.0 0.0 0.0; 0.0 16.0 0.0 42.0 20.0; 22.0 27.0 0.0 0.0 50.0; 7.0 4.0 0.0 0.0 19.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 -1.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0];
bg = [ -872.0; -739.0; -622.0; -533.0; -494.0; -596.0; -663.0; -891.0; -473.0; -273.0; -726.0; -774.0; -813.0; -1087.0; -809.0; -1369.0; -168.0; -975.0; -245.0; -98.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0];
bG = [];