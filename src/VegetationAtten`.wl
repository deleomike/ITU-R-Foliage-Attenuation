(* ::Package:: *)

BeginPackage["VegetationAtten`"]

Attenuation::usage := "Calculate specific attenuation of a tree."

Begin["`Private`"]

(*ASSUME: that thetai is the angle of entry for the transmitting signal,
and that thetas is the angle of exit for the receiver. thetai + thetas = pi
With respect to the vertical axis. Essentially theta is related 
to the heights of the transmitter and receiver.

ASSUME: phii is the azimuth for the entry, and phis is the azimuth of exit.
phii = phis because there is no change in direction assumed*)

Attenuation[A1_, A2_ , A3_ , F_ , Epsilon_, AngleEntry_, Azimuth_, BetaMax_, P_] :=

	Module[{},
		k = (2 * Pi * F) / (3 * 10 ^ 8);
(*		epsilon = 31 - I*8;*)
		epsilon = Epsilon;
		gamma = 0;
		
		thetas = AngleEntry;
		phis = Azimuth;
		
		phii = phis;
		thetai = Pi - thetas;
		
		ti=Sin[beta]*Sin[thetai]*Cos[alpha-phii]-Cos[beta]*Cos[thetai];
		
		lthetai = ArcCos[Cos[beta]*Cos[thetai]-Sin[beta]*Cos[alpha-phii]*Sin[thetai]];
		lthetas=lthetai/.thetai->Pi-thetas;
		lphii = ArcCos[(Sin[thetai]*(Cos[gamma]*Cos[beta]*Cos[alpha-phii]-Sin[gamma]*Sin[alpha-phii])+Cos[thetai]*Cos[gamma]*Sin[beta])/Sqrt[1-ti^2]];
		lphis=lphii/.thetai->(Pi-thetas);
		
		If[A1>A2 && A2>A3,
			(*Print["a"];*)
			e=Sqrt[1-(A2/A1)^2];
			KTemp=Integrate[1/Sqrt[1-e*Sin[x]^2],{x,0,Pi/2}];
			ETemp=Integrate[Sqrt[1-e*Sin[x]^2],{x,0,Pi/2}];
			g={
			(A3/A1)*Sqrt[1-e^2]*(KTemp-ETemp)/(e^2),
			(A3/A1)*(ETemp-(1-e^2)*KTemp)/(e^2*Sqrt[1-e^2]),
			1-(A3/A1)*ETemp/Sqrt[1-e^2]
			};
			qxalpha=k*(Sin[thetai]*Cos[phii-alpha]-Sin[thetas]*Cos[phis-alpha]);
			qyalpha=k*(Sin[thetai]*Sin[phii-alpha]-Sin[thetas]*Sin[phis-alpha]);
			qzalpha=k*(Cos[thetai]+Cos[thetas]);
			(*Qe=a1*Sqrt[(Cos[beta * qxalpha]+Sin[beta*qzalpha])^2+qyalpha^2]*)
			Qe = ((Cos[beta*qzalpha]+Sin[beta*qzalpha])^2*(A1^2*Cos[gamma]^2+A2^2*Sin[gamma]^2)+qyalpha*(Cos[beta*qzalpha]+Sin[beta*qzalpha])*Sin[2*gamma]*(A1^2-A2^2)+qyalpha^2*(A1*2*Sin[gamma]^2+A2^2*Cos[gamma]^2))^2;
			mu = 2*BesselJ[1,Qe]/Qe;
			
			att=1/((epsilon-1)*g+1);
			a33=Part[att,3];
			a22=Part[att,2];
			a11=Part[att,1];

			v0=(4*Pi/3)*A1*A2*A3;
			
			Fvv=-(k^2*v0/(4*Pi))*mu*(epsilon-1)*(a33*Sin[lthetai]*Sin[lthetas]-Cos[lthetai]*Cos[lthetas]*(a11*Cos[lphii]*Cos[lphis]+a22*Sin[lphii]*Sin[lphis]));
			Fhv=(k^2*v0/(4*Pi))*mu*(epsilon-1)*Cos[lthetai]*(a11*Sin[lphis]*Cos[lphii]-a22*Cos[lphis]*Sin[lphii]);
			Fvh=-(k^2*v0/(4*Pi))*mu*(epsilon-1)*Cos[lthetas]*(a11*Cos[lphis]*Sin[lphii]-a22*Sin[lphis]*Cos[lphii]);
			Fhh=(k^2*v0/(4*Pi))*mu*(epsilon-1)*(a11*Sin[lphis]*Sin[lphii]+a22*Cos[lphis]*Cos[lphii]);
			
			,
			If[A1==A2 && A1>A3,
				(*Print["b"];*)
				m=A1/A3;
				
				qxalpha=k*(Sin[thetai]*Cos[phii-alpha]-Sin[thetas]*Cos[phis-alpha]);
				qyalpha=k*(Sin[thetai]*Sin[phii-alpha]-Sin[thetas]*Sin[phis-alpha]);
				qzalpha=k*(Cos[thetai]+Cos[thetas]);
				Qe=A1*Sqrt[((Cos[beta*qxalpha]+Sin[beta*qzalpha])^2+qyalpha^2)];
				(*Qe = ((Cos[beta*qzalpha]+Sin[beta*qzalpha])^2*(A1^2*Cos[gamma]^2+A2^2*Sin[gamma]^2)+qyalpha*(Cos[beta*qzalpha]+Sin[beta*qzalpha])*Sin[2*gamma]*(A1^2-A2^2)+qyalpha^2*(A1*2*Sin[gamma]^2+A2^2*Cos[gamma]^2))^2;*)
				mu=2*BesselJ[1,Qe]/Qe;
				
				AT = ((1/(2*(m^2-1))) * ( ( m^2 / Sqrt[m^2-1] ) * ArcSin[(Sqrt[m^2-1]/m)]-1));
				AN = (m^2/(m^2-1))*(1-(1/Sqrt[m^2-1])*ArcSin[Sqrt[m^2-1]/m]);

				v0=(4*Pi/3)*A1*A2*A3;

				Fvv=(k^2 * (epsilon-1)/(4 * Pi)) * (v0) * (AN * Sin[lthetai] * Sin[lthetas] - AT * Cos[lthetai] * Cos[lthetas]*Cos[lphis-lphii])*mu;
				Fhv=(k^2 * (epsilon-1)/(4 * Pi)) * (v0) * Cos[lthetai]*Sin[lphis-lphii]*AT*mu;
				Fvh=(k^2 * (epsilon-1)/(4 * Pi)) * (v0) * Cos[lthetas]*Sin[lphis-lphii]*AT*mu;
				Fhh=(k^2 * (epsilon-1)/(4 * Pi)) * (v0) * Cos[lphis-lphii]*AT*mu;
				,
				If[A1==A2 && A1<A3,
					(*Print["c"];*)
					
					b = Sqrt[1 - (A1/A3)^2];
					
					g1 = (b * (b^2 -  1) / 2) * ( (b/(b^2-1)) + (1/2) * Log[(b-1)/(b+1)]);
					g3 = -(b^2 - 1) * ( (b/2) * Log[(b - 1)/(b + 1)] + 1);
					
					AT=1/((epsilon-1)*g1+1);
					AN=1/((epsilon-1)*g3+1);
					
					qxalpha=k*(Sin[thetai]*Cos[phii-alpha]-Sin[thetas]*Cos[phis-alpha]);
					qzalpha=k*(Cos[thetai]+Cos[thetas]);
					qz = Sin[beta * qxalpha] - Cos[beta * qzalpha];
					mu = Sin[qz * A3] / (qz * A3);
					
					Fvv=(k^2 * (epsilon-1)/(4 * Pi)) * (v0) * (AN * Sin[lthetai] * Sin[lthetas] - AT * Cos[lthetai] * Cos[lthetas]*Cos[lphis-lphii])*mu;
					Fhv=(k^2 * (epsilon-1)/(4 * Pi)) * (v0) * Cos[lthetai]*Sin[lphis-lphii]*AT*mu;
					Fvh=(k^2 * (epsilon-1)/(4 * Pi)) * (v0) * Cos[lthetas]*Sin[lphis-lphii]*AT*mu;
					Fhh=(k^2 * (epsilon-1)/(4 * Pi)) * (v0) * Cos[lphis-lphii]*AT*mu;
					
					,
					0
					];
				];
			];
			
		tvi =-(Sin[beta]*Cos[thetai]*Cos[alpha-phii]+Cos[beta]*Sin[thetai]);
		thi = Sin[beta]*Sin[alpha-phii];
		tvs = -(Sin[beta]*Cos[Pi-thetas]*Cos[alpha-phis]+Cos[beta]*Sin[Pi-thetas]);
		ths = Sin[beta]*Sin[alpha-phis];
		d = Sqrt[(tvs^2+ths^2)*(tvi^2+thi^2)];
		FVV = (1/d)*(tvs*(Fvv*tvi-Fvh*thi)-ths*(Fhv*tvi+Fhh*thi));
	(*	Print[FVV];*)
	(*	Print["--"];*)
	(*	Print[FVV];*)
(*		FVV = FVV/.thetai->AngleEnter/.phii->(AngleExit);*)
(*		Print[AngleEnter];*)
(*		Print[FVV];*)
(*		Print[FVV/.beta->(AngleEnter)/.alpha->AngleExit];*)
		
(*		Print[Fvv];
		Print[FVV];*)
		EField = NIntegrate[(Sin[beta]/(Cos[0]-Cos[BetaMax]))* FVV,{beta,0,BetaMax},{alpha,0,2*Pi}];
		Feq= P * EField;
		KP = k +(2*Pi/k) * Feq;
		alphac=-8.686*Im[KP]
		]
	
End[]

EndPackage[]

































