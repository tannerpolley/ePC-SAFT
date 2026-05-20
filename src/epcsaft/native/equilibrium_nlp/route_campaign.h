#pragma once

#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft::native::equilibrium_nlp {

template <typename Result, typename Attempt>
class RouteCampaign {
public:
    using QualityFunction = std::function<int(const Result&)>;

    RouteCampaign(
        QualityFunction quality,
        std::string initial_point_strategy,
        std::string empty_campaign_seed_name
    )
        : quality_(std::move(quality)),
          initial_point_strategy_(std::move(initial_point_strategy)),
          empty_campaign_seed_name_(std::move(empty_campaign_seed_name)) {}

    void reserve(std::size_t count) {
        attempts_.reserve(count);
    }

    void record_attempt(const Result& result, Attempt attempt) {
        attempts_.push_back(std::move(attempt));
        if (!have_best_ || quality_(result) > quality_(best_)) {
            best_ = result;
            have_best_ = true;
        }
    }

    Result finish_result() const {
        Result out = best_;
        if (!have_best_) {
            out.initial_point_strategy = initial_point_strategy_;
            out.seed_name = empty_campaign_seed_name_;
        }
        out.seed_attempts = attempts_;
        return out;
    }

private:
    QualityFunction quality_;
    std::string initial_point_strategy_;
    std::string empty_campaign_seed_name_;
    Result best_;
    bool have_best_ = false;
    std::vector<Attempt> attempts_;
};

}  // namespace epcsaft::native::equilibrium_nlp
