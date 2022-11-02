Fx = [ 49.0; -14.0; 3.0; -20.0; -23.0; 17.0; -26.0; 22.0; 4.0; -28.0; -24.0; 44.0; -14.0; -7.0; -11.0];
Fy = [ -14.0; -19.0; 47.0; 17.0; 42.0];
Gx = [];
Gy = [];
gx = [ 0.0 46.0 0.0 43.0 35.0 19.0 0.0 0.0 42.0 44.0 40.0 31.0 0.0 10.0 0.0; 50.0 4.0 0.0 0.0 44.0 26.0 0.0 47.0 0.0 24.0 40.0 34.0 8.0 0.0 0.0; 41.0 32.0 13.0 17.0 8.0 29.0 26.0 0.0 0.0 23.0 0.0 0.0 38.0 0.0 1.0; 8.0 31.0 9.0 16.0 4.0 5.0 32.0 48.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 48.0 0.0 6.0 33.0 21.0 18.0 0.0 49.0 0.0 0.0 2.0 5.0 3.0 0.0; 0.0 0.0 15.0 14.0 46.0 12.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 22.0 0.0; 0.0 18.0 46.0 2.0 9.0 0.0 0.0 0.0 0.0 43.0 21.0 28.0 0.0 40.0 39.0; 45.0 34.0 22.0 6.0 26.0 28.0 33.0 8.0 25.0 40.0 30.0 22.0 28.0 7.0 16.0; 46.0 10.0 48.0 45.0 3.0 0.0 35.0 45.0 11.0 0.0 4.0 0.0 46.0 25.0 0.0; 5.0 34.0 1.0 44.0 14.0 1.0 7.0 0.0 43.0 27.0 45.0 1.0 8.0 26.0 5.0; 0.0 7.0 17.0 30.0 0.0 1.0 21.0 36.0 0.0 0.0 43.0 48.0 6.0 3.0 0.0; 47.0 15.0 0.0 3.0 30.0 10.0 14.0 45.0 21.0 0.0 22.0 0.0 40.0 47.0 0.0; 0.0 0.0 1.0 45.0 24.0 0.0 25.0 0.0 36.0 0.0 0.0 0.0 25.0 0.0 0.0; 2.0 29.0 19.0 0.0 2.0 0.0 5.0 0.0 46.0 0.0 30.0 0.0 28.0 31.0 0.0; 13.0 10.0 0.0 7.0 7.0 49.0 10.0 30.0 38.0 44.0 31.0 29.0 11.0 12.0 15.0; 0.0 37.0 7.0 0.0 1.0 0.0 0.0 40.0 45.0 37.0 36.0 39.0 38.0 0.0 0.0; 9.0 6.0 0.0 24.0 0.0 38.0 40.0 20.0 42.0 36.0 44.0 14.0 5.0 0.0 26.0; 10.0 40.0 15.0 8.0 49.0 30.0 26.0 48.0 0.0 50.0 44.0 42.0 0.0 43.0 49.0; 17.0 14.0 45.0 25.0 0.0 26.0 25.0 42.0 46.0 6.0 5.0 35.0 28.0 29.0 31.0; 0.0 26.0 0.0 0.0 9.0 18.0 24.0 25.0 46.0 19.0 22.0 30.0 21.0 30.0 10.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0];
gy = [ 5.0 1.0 0.0 42.0 21.0; 0.0 7.0 0.0 0.0 21.0; 6.0 24.0 0.0 4.0 34.0; 28.0 45.0 36.0 0.0 35.0; 40.0 0.0 3.0 30.0 19.0; 35.0 2.0 0.0 0.0 0.0; 32.0 47.0 0.0 5.0 0.0; 0.0 19.0 46.0 0.0 3.0; 28.0 32.0 28.0 22.0 23.0; 3.0 41.0 17.0 11.0 42.0; 46.0 32.0 17.0 45.0 0.0; 9.0 7.0 36.0 50.0 11.0; 0.0 29.0 41.0 27.0 0.0; 0.0 31.0 10.0 30.0 0.0; 29.0 7.0 13.0 48.0 32.0; 7.0 50.0 11.0 0.0 44.0; 8.0 28.0 7.0 9.0 34.0; 45.0 46.0 44.0 36.0 35.0; 11.0 13.0 29.0 13.0 47.0; 16.0 3.0 21.0 50.0 26.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 -1.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0];
bg = [ -501.0; -609.0; -513.0; -416.0; -218.0; -78.0; -225.0; -592.0; -576.0; -407.0; -484.0; -381.0; -235.0; -56.0; -638.0; -478.0; -687.0; -761.0; -686.0; -458.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0];
bG = [];
