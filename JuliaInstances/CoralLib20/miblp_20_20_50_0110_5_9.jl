Fx = [ 20.0; 18.0; 32.0; -14.0; -1.0; 11.0; 22.0; -29.0; 24.0; 33.0; -41.0; 12.0; -30.0; -9.0; -6.0];
Fy = [ 35.0; 11.0; -29.0; 3.0; -41.0];
Gx = [];
Gy = [];
gx = [ 36.0 19.0 0.0 17.0 11.0 0.0 33.0 8.0 4.0 0.0 47.0 50.0 6.0 32.0 33.0; 38.0 12.0 21.0 12.0 47.0 19.0 46.0 7.0 0.0 0.0 9.0 38.0 8.0 39.0 19.0; 7.0 28.0 0.0 21.0 14.0 49.0 0.0 41.0 0.0 0.0 26.0 0.0 4.0 12.0 22.0; 35.0 0.0 10.0 20.0 12.0 50.0 6.0 25.0 13.0 33.0 18.0 17.0 23.0 2.0 47.0; 6.0 0.0 48.0 0.0 41.0 0.0 43.0 1.0 6.0 27.0 0.0 41.0 37.0 0.0 30.0; 0.0 0.0 41.0 41.0 0.0 3.0 0.0 0.0 0.0 0.0 31.0 41.0 11.0 50.0 45.0; 44.0 11.0 37.0 45.0 26.0 0.0 8.0 15.0 2.0 18.0 0.0 6.0 49.0 27.0 21.0; 2.0 0.0 3.0 27.0 0.0 30.0 33.0 19.0 49.0 36.0 5.0 2.0 0.0 0.0 14.0; 0.0 44.0 0.0 21.0 29.0 0.0 0.0 0.0 0.0 50.0 0.0 0.0 0.0 36.0 3.0; 0.0 0.0 28.0 0.0 24.0 14.0 50.0 0.0 0.0 3.0 0.0 35.0 36.0 17.0 0.0; 40.0 16.0 38.0 40.0 50.0 47.0 30.0 27.0 15.0 28.0 0.0 0.0 9.0 9.0 25.0; 17.0 8.0 14.0 7.0 8.0 47.0 8.0 43.0 10.0 0.0 0.0 1.0 32.0 14.0 19.0; 37.0 36.0 12.0 0.0 41.0 0.0 31.0 32.0 2.0 0.0 34.0 0.0 0.0 30.0 40.0; 8.0 19.0 46.0 17.0 23.0 11.0 11.0 19.0 16.0 13.0 19.0 45.0 1.0 19.0 21.0; 36.0 37.0 1.0 0.0 0.0 19.0 39.0 18.0 5.0 26.0 0.0 34.0 3.0 42.0 27.0; 36.0 16.0 19.0 16.0 4.0 47.0 28.0 19.0 8.0 45.0 0.0 22.0 8.0 13.0 33.0; 0.0 0.0 48.0 13.0 35.0 0.0 47.0 0.0 28.0 31.0 25.0 0.0 31.0 29.0 2.0; 0.0 46.0 0.0 0.0 8.0 0.0 0.0 0.0 27.0 48.0 0.0 31.0 30.0 15.0 8.0; 32.0 0.0 2.0 23.0 29.0 27.0 50.0 42.0 35.0 5.0 19.0 0.0 22.0 44.0 42.0; 3.0 32.0 32.0 48.0 0.0 26.0 0.0 39.0 24.0 31.0 11.0 9.0 23.0 12.0 33.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0];
gy = [ 0.0 29.0 12.0 0.0 49.0; 40.0 18.0 49.0 13.0 30.0; 7.0 28.0 36.0 48.0 17.0; 21.0 0.0 11.0 30.0 24.0; 2.0 49.0 18.0 14.0 0.0; 27.0 0.0 0.0 4.0 20.0; 22.0 0.0 0.0 15.0 28.0; 0.0 0.0 0.0 50.0 0.0; 0.0 0.0 0.0 0.0 32.0; 0.0 20.0 0.0 13.0 21.0; 4.0 32.0 50.0 0.0 36.0; 12.0 40.0 9.0 25.0 4.0; 0.0 0.0 35.0 40.0 20.0; 12.0 21.0 39.0 27.0 2.0; 0.0 0.0 0.0 21.0 28.0; 24.0 22.0 41.0 0.0 22.0; 17.0 41.0 34.0 31.0 22.0; 0.0 0.0 23.0 0.0 12.0; 18.0 12.0 41.0 13.0 4.0; 13.0 0.0 0.0 0.0 39.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 -1.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0];
bg = [ -724.0; -934.0; -456.0; -566.0; -609.0; -907.0; -706.0; -113.0; -510.0; -436.0; -737.0; -360.0; -669.0; -788.0; -672.0; -769.0; -684.0; -564.0; -592.0; -757.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0];
bG = [];