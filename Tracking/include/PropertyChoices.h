/*
 * Copyright (c) 2020-2026 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef PROPERTY_CHOICES_H
#define PROPERTY_CHOICES_H

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

/** @class PropertyChoices
 *
 *  Minimal stand-in for the `choices` argument of Python's argparse.add_argument(), for
 *  string-valued Gaudi properties that may only take one of a fixed set of values. It replaces a
 *  hand-written if/else chain by a single table that is also the one source of truth for the list
 *  of valid values quoted in the property documentation and in error messages.
 *
 *  Gaudi itself has no equivalent. Gaudi::Property does take a VERIFIER template parameter, but
 *  the only two verifiers it ships are Gaudi::Details::Property::NullVerifier and
 *  Gaudi::Details::Property::BoundedVerifier (numeric lower/upper bounds, exposed as
 *  Gaudi::CheckedProperty). A verifier is moreover default-constructed by the property and already
 *  invoked on the default value inside the property constructor, so there is no clean way to teach
 *  one a list of allowed strings from the owning algorithm.
 *
 *  Usage: declare the choices next to the property, then translate the configured string once in
 *  initialize() and keep the enum for use in the event loop.
 *
 *      enum class Mode { Fast, Slow };
 *      static constexpr auto s_modeChoices =
 *          makePropertyChoices<Mode>(std::pair{"Fast", Mode::Fast}, std::pair{"Slow", Mode::Slow});
 *
 *      Gaudi::Property<std::string> m_mode{this, "Mode", "Fast",
 *                                          "Which mode to run, one of " + s_modeChoices.list()};
 *
 *      StatusCode initialize() override {
 *        const auto mode = s_modeChoices.parse(m_mode);
 *        if (!mode) {
 *          error() << "Invalid Mode '" << m_mode.value() << "', expected one of "
 *                  << s_modeChoices.list() << "." << endmsg;
 *          return StatusCode::FAILURE;
 *        }
 *        m_parsedMode = *mode;
 *        return StatusCode::SUCCESS;
 *      }
 *
 *  @author Andreas Loeschcke Centeno
 */
template <typename ENUM, std::size_t N>
class PropertyChoices {
public:
  using Choice = std::pair<std::string_view, ENUM>;

  constexpr explicit PropertyChoices(std::array<Choice, N> choices) : m_choices(choices) {}

  /// Translate one of the allowed strings into its enum value, or return std::nullopt if the
  /// string is not one of the choices
  constexpr std::optional<ENUM> parse(std::string_view value) const {
    for (const auto& choice : m_choices) {
      if (choice.first == value)
        return choice.second;
    }
    return std::nullopt;
  }

  /// The allowed values rendered as `'A', 'B', 'C'`, for property documentation and error messages
  std::string list() const {
    std::string rendered;
    for (const auto& choice : m_choices) {
      if (!rendered.empty())
        rendered += ", ";
      rendered += '\'';
      rendered += choice.first;
      rendered += '\'';
    }
    return rendered;
  }

private:
  std::array<Choice, N> m_choices;
};

/// Build a PropertyChoices out of `std::pair{"Name", Enum::Value}` entries, deducing the number of
/// choices so that it never has to be kept in sync by hand
template <typename ENUM, typename... NAMES>
constexpr auto makePropertyChoices(std::pair<NAMES, ENUM>... choices) {
  return PropertyChoices<ENUM, sizeof...(choices)>{std::array<std::pair<std::string_view, ENUM>, sizeof...(choices)>{
      std::pair<std::string_view, ENUM>{choices.first, choices.second}...}};
}

#endif // PROPERTY_CHOICES_H
