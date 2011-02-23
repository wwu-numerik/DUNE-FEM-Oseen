clear all;
run_matlab;

p_dofs_count = size(p_exakt,1);
u_dofs_count = size(u_exakt,1);

% allg matrizen
A = Y - X * M_invers * W;
E = E .* -1;


% gewoehnliches schurkomplement
A_inv = inv(A);
F = H2 - X * M_invers * H1;
schur = R + E * A_inv * Z ;
schur_F = E * A_inv * F - H3;
p_s = schur \ schur_F;
u_s = ( A \ ( F - Z * p_s ) ) ;

% alt. schurkomplement
%  C_inv = inv(-R);
%  schur_alt = ( A - Z * C_inv * E );
%  schur_F_alt = F + ( Z * C_inv * H3 );
%  o_u = schur_alt / schur_F_alt';
%  o_p = C_inv * ( H3 + E * o_u);

% vollsystem
S = [ A Z ; E R ] ;
SF = [ F ; H3 ];
up = S \ SF;
up_u = up(1:u_dofs_count);
up_p = up(u_dofs_count:u_dofs_count+p_dofs_count-1);

%up_exakt = [ u_exakt ; p_exakt ];
%diff = up - up_exakt;
%diff_s = up_s - up_exakt;
%norm ( diff_s )
%p_diff = p_exakt - up(size(u_exakt)+1:size(up) );
%u_diff = u_exakt - up(1:size(u_exakt));
%p_diff_norm = norm( p_diff )
%u_diff_norm = norm( u_diff )
%up_diff_norm = norm( diff )
%  rank_S = sprank( S )
%  rank_A = sprank( A )
%  rank_Y = sprank( Y )
%  rank_E_neg = sprank( E_neg )
%  rank_C = sprank( C )
%  rank_Z = sprank( Z )
%  cond_S = condest( S )
%  cond_A = condest( A )
%  cond_Y = condest( Y )
%  cond_C = condest( C );
%  chol ( A );
%  size( A )[0] == rank_A
%cond_Z = condest( Z )
%  cond_E_neg = condest( E_neg )

p_dof_ids = [ 1 : p_dofs_count ];
u_dof_ids = [ 1 : u_dofs_count ];
figure(1);
%  plot( p_dof_ids, p_computed, 'go', p_dof_ids, up_p, 'rx');
%  p_s = p_s + (ones( p_dofs_count, 1 ) .* 15.6655);
plot( p_dof_ids, p_computed, 'go', p_dof_ids, p_s, 'rx');
%  plot( i, up_p, 'go');
%  plot( i, DD, 'go', i, p_computed, 'rx');
xlabel('dof');
ylabel('value');
title('pressure');

figure(2);
plot( u_dof_ids, u_computed, 'bo', u_dof_ids, up_u, 'rx');
%  plot( u_dof_ids, u_computed, 'bo', u_dof_ids, u_s, 'rx');
%  plot( i, p_s, 'rx');
xlabel('dof');
ylabel('value');
title('velocity')
