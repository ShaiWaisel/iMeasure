a = [0 -2 0;1 -2 0; 1 -2 1; 0 -2 0];
b = [10 0 0; 10 1 0; 10 1.1 -1; 10 0 0];
v = [0.5 0 0.5 1];
[r,f,e]=absor(a',b')
q=v*r.M'
figure, 
plot3(a(:,1), a(:,2), a(:,3),'r'); 
hold on;
axis equal;
plot3(b(:,1), b(:,2), b(:,3),'b'); 
scatter3(v(1),v(2),v(3),'red');
scatter3(q(1),q(2),q(3),'ob','filled')
for i=1:size(a,1)
    text(a(i,1),a(i,2),a(i,3)+0.1,num2str(i),'Color','red');
    text(b(i,1),b(i,2),b(i,3)+0.1,num2str(i),'Color','blue');
end

