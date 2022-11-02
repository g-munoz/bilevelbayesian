Fx = [ -40.0; -18.0; 22.0; 18.0; 4.0];
Fy = [ 36.0; 36.0; -22.0; 3.0; 1.0; 25.0; 24.0; 31.0; -30.0; 39.0; -30.0; -3.0; -35.0; -39.0; 14.0];
Gx = [];
Gy = [];
gx = [ 4.0 48.0 0.0 38.0 1.0; 31.0 0.0 8.0 0.0 0.0; 0.0 6.0 44.0 40.0 0.0; 0.0 31.0 0.0 0.0 0.0; 27.0 50.0 24.0 0.0 5.0; 8.0 19.0 19.0 0.0 40.0; 38.0 0.0 6.0 47.0 0.0; 18.0 44.0 43.0 10.0 35.0; 0.0 39.0 11.0 16.0 30.0; 9.0 24.0 11.0 11.0 30.0; 0.0 8.0 42.0 0.0 33.0; 2.0 0.0 22.0 0.0 0.0; 50.0 0.0 31.0 35.0 38.0; 0.0 41.0 0.0 0.0 0.0; 0.0 43.0 10.0 10.0 22.0; 28.0 32.0 3.0 45.0 48.0; 34.0 40.0 31.0 33.0 0.0; 0.0 0.0 13.0 0.0 0.0; 0.0 19.0 0.0 19.0 19.0; 28.0 12.0 0.0 9.0 17.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 -1.0];
gy = [ 0.0 4.0 13.0 23.0 43.0 28.0 33.0 0.0 41.0 46.0 15.0 47.0 5.0 25.0 36.0; 26.0 42.0 0.0 0.0 0.0 11.0 0.0 0.0 5.0 4.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 45.0 0.0 47.0 27.0 0.0 22.0 30.0 0.0 42.0 36.0 9.0 0.0 21.0; 0.0 29.0 10.0 0.0 34.0 25.0 38.0 0.0 14.0 31.0 3.0 0.0 23.0 11.0 0.0; 50.0 27.0 7.0 16.0 37.0 49.0 39.0 26.0 10.0 40.0 19.0 31.0 21.0 42.0 2.0; 35.0 5.0 0.0 44.0 0.0 0.0 0.0 0.0 0.0 36.0 38.0 0.0 19.0 0.0 0.0; 0.0 29.0 27.0 18.0 37.0 13.0 6.0 0.0 0.0 1.0 24.0 0.0 8.0 28.0 0.0; 25.0 25.0 19.0 3.0 47.0 0.0 47.0 37.0 31.0 22.0 17.0 6.0 8.0 46.0 41.0; 23.0 26.0 25.0 2.0 33.0 43.0 0.0 48.0 9.0 0.0 35.0 14.0 5.0 48.0 1.0; 26.0 24.0 11.0 0.0 35.0 11.0 0.0 21.0 11.0 21.0 41.0 8.0 2.0 37.0 49.0; 0.0 0.0 0.0 35.0 0.0 34.0 20.0 0.0 0.0 14.0 0.0 0.0 0.0 25.0 19.0; 0.0 26.0 8.0 45.0 46.0 23.0 0.0 38.0 49.0 0.0 0.0 20.0 22.0 38.0 0.0; 46.0 33.0 17.0 17.0 2.0 0.0 0.0 36.0 0.0 37.0 5.0 44.0 22.0 46.0 32.0; 23.0 0.0 6.0 5.0 27.0 33.0 19.0 40.0 36.0 34.0 8.0 0.0 43.0 19.0 23.0; 0.0 37.0 47.0 0.0 47.0 0.0 0.0 9.0 19.0 1.0 0.0 35.0 5.0 0.0 0.0; 19.0 31.0 50.0 24.0 44.0 35.0 0.0 40.0 0.0 0.0 2.0 42.0 7.0 4.0 49.0; 40.0 13.0 38.0 9.0 13.0 0.0 17.0 0.0 48.0 26.0 24.0 0.0 11.0 0.0 25.0; 0.0 0.0 26.0 0.0 22.0 30.0 0.0 27.0 0.0 12.0 25.0 0.0 0.0 38.0 25.0; 27.0 18.0 0.0 0.0 0.0 13.0 6.0 0.0 45.0 25.0 0.0 13.0 0.0 0.0 45.0; 34.0 23.0 0.0 0.0 46.0 0.0 49.0 0.0 15.0 48.0 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
bg = [ -485.0; -63.0; -439.0; -180.0; -368.0; -187.0; -134.0; -522.0; -329.0; -511.0; -203.0; -231.0; -262.0; -422.0; -66.0; -382.0; -439.0; -312.0; -425.0; -177.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0];
bG = [];
