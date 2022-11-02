Fx = [ 6.0; -47.0; 26.0; 27.0; 27.0];
Fy = [ 8.0; -44.0; 14.0; -15.0; 19.0; -3.0; 14.0; -26.0; 31.0; 26.0];
Gx = [];
Gy = [];
gx = [ 0.0 43.0 48.0 10.0 36.0; 0.0 0.0 0.0 11.0 4.0; 2.0 0.0 2.0 25.0 39.0; 47.0 0.0 0.0 12.0 29.0; 4.0 50.0 45.0 0.0 25.0; 0.0 0.0 0.0 0.0 14.0; 0.0 5.0 43.0 38.0 0.0; 11.0 35.0 28.0 0.0 17.0; 34.0 30.0 47.0 3.0 41.0; 16.0 37.0 0.0 12.0 42.0; 49.0 44.0 32.0 7.0 0.0; 14.0 23.0 13.0 49.0 20.0; 0.0 0.0 25.0 0.0 9.0; 0.0 31.0 11.0 27.0 0.0; 41.0 28.0 0.0 0.0 42.0; 0.0 9.0 16.0 0.0 44.0; 11.0 7.0 0.0 28.0 14.0; 29.0 42.0 26.0 24.0 18.0; 34.0 26.0 38.0 21.0 40.0; 15.0 8.0 0.0 26.0 44.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 -1.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0];
gy = [ 10.0 1.0 0.0 0.0 0.0 27.0 0.0 0.0 0.0 47.0; 0.0 0.0 0.0 0.0 41.0 28.0 0.0 0.0 41.0 19.0; 0.0 0.0 0.0 0.0 0.0 44.0 0.0 0.0 0.0 42.0; 4.0 16.0 36.0 5.0 41.0 25.0 49.0 35.0 41.0 22.0; 39.0 0.0 13.0 0.0 0.0 40.0 0.0 0.0 45.0 42.0; 1.0 0.0 15.0 7.0 26.0 0.0 0.0 45.0 44.0 8.0; 0.0 0.0 25.0 0.0 32.0 0.0 0.0 11.0 0.0 45.0; 5.0 10.0 28.0 7.0 16.0 0.0 33.0 26.0 15.0 5.0; 8.0 25.0 30.0 0.0 21.0 47.0 25.0 25.0 0.0 1.0; 2.0 16.0 5.0 0.0 50.0 0.0 41.0 15.0 40.0 37.0; 0.0 44.0 0.0 11.0 16.0 8.0 0.0 21.0 18.0 5.0; 19.0 43.0 23.0 6.0 23.0 0.0 0.0 4.0 50.0 28.0; 0.0 45.0 39.0 2.0 4.0 5.0 21.0 44.0 0.0 14.0; 3.0 34.0 22.0 21.0 24.0 9.0 24.0 6.0 40.0 38.0; 42.0 17.0 0.0 34.0 0.0 13.0 8.0 44.0 0.0 16.0; 31.0 28.0 24.0 0.0 0.0 0.0 0.0 27.0 44.0 0.0; 25.0 0.0 0.0 0.0 20.0 31.0 33.0 0.0 0.0 35.0; 12.0 6.0 22.0 50.0 48.0 8.0 40.0 23.0 42.0 43.0; 38.0 39.0 0.0 21.0 25.0 34.0 13.0 23.0 36.0 44.0; 40.0 0.0 0.0 48.0 49.0 0.0 18.0 38.0 40.0 27.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0];
bg = [ -187.0; -93.0; -209.0; -293.0; -134.0; -196.0; -329.0; -182.0; -172.0; -272.0; -173.0; -340.0; -225.0; -290.0; -244.0; -117.0; -250.0; -380.0; -324.0; -362.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0; -1500.0; 0.0];
bG = [];
