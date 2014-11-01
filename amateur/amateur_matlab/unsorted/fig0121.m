%compare Lprime and L from ahfout
%ahfinput

%lambda,lambdaE describe \lambda' (Bullock 2001) and \lambda (Peebles)

summary(dataset(lambda))
summary(dataset(lambdaE))
sum(isnan(lambdaE))
rose(log10(lambdaE),100)
hold on;
rose(log10(lambda),100)
hold off;