function plotting()
x = [-1:0.01:1];
for i = 1:length( x )
    for j = 1:length( x )
        x_1 = x( i );
        x_2 = x( j );
        p_exact( i, j ) = exp( x_1 ) * x_2 * sin( x_2 );
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% plot the results                                                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
screen_size = get(0,'ScreenSize');
figure( 'Position', [1 1 screen_size(3) screen_size(4)],...
        'Name', 'exact pressure' );
mesh( x, x, p_exact );
title('exact pressure');
xlabel('-1 < x < 1');
ylabel('-1 < y < 1');
legend('exact pressure');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %