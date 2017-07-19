void scaling_(Par *par, Ellipsoid *ellipsoid1, Ellipsoid *ellipsoid2, SizeEllipsoid *sizes, Overlap *overlap){
	// Updates the contact points between ellipsoid1 and ellipsoid2.

	// DETERMINATION OF THE INVERSE OF THE REDUCED BELONGING MATRIX

	double Bainv[9]; // inverse of the reduced belong matrix of ellipsoid1
	inv_red_bel(&Bainv[0],ellipsoid1,sizes);
	printf("\nBainv:\n");
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			printf("%f ",Bainv[i*3+j]);
		}
		printf("\n");
	}

	double Bbinv[9]; // inverse of the reduced belong matrix of ellipsoid2
	inv_red_bel(&Bbinv[0],ellipsoid2,sizes);
	printf("\nBbinv:\n");
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			printf("%f ",Bbinv[i*3+j]);
		}
		printf("\n");
	}

	// DETERMINATION OF POLYNOMIAL Y_{AB}

	Pol yab[9]; // polynomial Y_{AB}
	y_ab(&yab[0],&Bainv[0],&Bbinv[0]);
	printf("\nyab:\n");
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			for (int c = 0; c <= yab[i*3+j].deg; c++){
				printf("%f ",yab[i*3+j].coef[c]);
			}
			printf(" || ");
		}
		printf("\n");
	}

	// DETERMINATION OF THE ADJOINT AND DETERMINANT OF Y_{AB}

	Pol adj[9]; // adjoint of matrix Y_{AB}
	adj_mat3_pol(&adj[0],yab);
	printf("\nadj:\n");
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			for (int c = 0; c <= adj[i*3+j].deg; c++){
				printf("%f ",adj[i*3+j].coef[c]);
			}
			printf(" || ");
		}
		printf("\n");
	}

	Pol qab; // determinant of matrix Y_{AB} (polynomial Q_{AB})
	det3_mat_pol(&qab,&yab[0]);
	printf("\nqab:\n");
	for (int c = 0; c <= qab.deg; c++){
		printf("%f ",qab.coef[c]);
	}
	printf("\n");

	// DETERMINATION OF R_{AB}

	double rab[3]; // \vec{r_B} - \vec{r_A}
	matrix_diff(&rab[0],ellipsoid2->r,ellipsoid1->r,3);
	printf("\nrab:\n");
	for (int c = 0; c < 3; c++){
		printf("%f ",rab[c]);
	}
	printf("\n");

	// DETERMINATION OF POLYNOMIAL P_{AB}

	Pol pab; // polynomial P_{AB}
	p_ab(&pab,&adj[0],&rab[0]);
	printf("\npab:\n");
	for (int c = 0; c <= pab.deg; c++){
		printf("%f ",pab.coef[c]);
	}
	printf("\n");

	// DETERMINATION OF POLYNOMIAL H_{AB}

	Pol hab; // polynomial H_{AB}
	hab.deg = 7; // H_{AB} is a polynomial of degree 7
	hab.coef = malloc(sizeof(double) * 8);
	init(&hab.coef[0],8); // initialisation of the coefficients of the polynomial

	for (int copp = 1; copp <= pab.deg; copp++){
		for (int coq = 0; coq <= qab.deg; coq++){
			hab.coef[copp + coq - 1] += copp*pab.coef[copp]*qab.coef[coq]; // P'_{AB} Q_{AB}
		}
	}

	for (int cop = 0; cop <= pab.deg; cop++){
		for (int coqq = 1; coqq <= qab.deg; coqq++){
			hab.coef[cop + coqq - 1] -= pab.coef[cop]*coqq*qab.coef[coqq]; // - P_{AB} Q'_{AB}
		}
	}

	printf("\nhab:\n");
	for (int c = 0; c <= hab.deg; c++){
		printf("%f ",hab.coef[c]);
	}
	printf("\n");

	// DETERMINATION OF THE ROOT OF H_{AB} IN [0,1]

	double root = (sizes + ellipsoid1->size)->R/((sizes + ellipsoid1->size)->R + (sizes + ellipsoid2->size)->R); // root of the polynomial in [0,1]
	root = root_Newton_pol(&hab,root,par->epsilon);
	printf("\nroot: %f\n",root);

	// DETERMINATION OF THE RESCALING FACTOR

	overlap->mu = sqrt(eval_pol(&pab,root)/eval_pol(&qab,root)); // mu = sqrt(potential)
	printf("\npab:\n");
	for (int c = 0; c <= pab.deg; c++){
		printf("%f ",pab.coef[c]);
	}
	printf("\n");
	printf("pab(root) = %f\n",eval_pol(&pab,root));
	printf("\nqab:\n");
	for (int c = 0; c <= qab.deg; c++){
		printf("%f ",qab.coef[c]);
	}
	printf("\n");
	printf("qab(root) = %f\n",eval_pol(&qab,root));

	printf("\nsqrt(pab(root)/qab(root)) = %f\n",sqrt(eval_pol(&pab,root)/eval_pol(&qab,root)));


	// CONTACT POINTS

	if (overlap->mu < 1){
		// ellipsoids are overlapping
		overlap->overlapping = 1;

		double det = eval_pol(&qab,root); // determinant of Y_{AB}

		double yabinvrab[3]; // Y_{AB}^-1 (\vec{r_B} - \vec{r_A})
		init(&yabinvrab[0],3); // initialisation
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				yabinvrab[i] += eval_pol(&adj[i * 3 + j],root) * rab[j] / det; // adj(Y_{AB}) (\vec{r_B} - \vec{r_A}) / det(Y_{AB}) = Y_{AB}^-1 (\vec{r_B} - \vec{r_A})
			}
		}

		for (int i = 0; i < 3; i++){
			// contact point
			overlap->r_c[i] = ellipsoid2->r[i]; // r_B
			for (int j = 0; j < 3; j++){
				overlap->r_c[i] -= root * Bbinv[i * 3 + j] * yabinvrab[j]; // - lambda B_B^-1 Y_{AB}^-1 (\vec{r_B} - \vec{r_A})
			}
			printf("%f\n",overlap->r_c[i]);
			// contact point on particles
			overlap->r_c1[i] = ellipsoid1->r[i] + (overlap->r_c[i] - ellipsoid1->r[i])/overlap->mu; // contact point on particle 1
			overlap->r_c2[i] = ellipsoid2->r[i] + (overlap->r_c[i] - ellipsoid2->r[i])/overlap->mu; // contact point on particle 2
		}
	}
	else {
		// particles are non-overlapping
		overlap->overlapping = 0;
	}

	// FREE POLYNOMIALS

	free_mat_pol(&yab[0],9);
	free_mat_pol(&adj[0],9);
	free(qab.coef);
	free(pab.coef);
	free(hab.coef);

}
