using ShapML

ShapML.shap(model::MLJBase.Machine, X=first(mach.args).data;
            target_features=nothing, sample_size::Int=60,
            reference=nothing, parallel::Symbol=:features) =
    ShapML.shap(; explain=X, reference, model, target_features, sample_size,
                parallel, predict_function=
                    (m, x) -> DataFrame(y=predict(m, x)))
