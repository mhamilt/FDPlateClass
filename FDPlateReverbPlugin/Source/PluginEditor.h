/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class FdplateReverbPluginAudioProcessorEditor  : public AudioProcessorEditor, private Slider::Listener
{
public:
    FdplateReverbPluginAudioProcessorEditor (FdplateReverbPluginAudioProcessor&);
    ~FdplateReverbPluginAudioProcessorEditor();

    //==============================================================================
    void paint (Graphics&) override;
    void resized() override;
    void sliderValueChanged (Slider* slider) override;
private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    Slider gainSlider;
    FdplateReverbPluginAudioProcessor& processor;
    ScopedPointer<AudioProcessorValueTreeState::SliderAttachment> gainSlideAttach;
//    std::unique_ptr<AudioProcessorValueTreeState::SliderAttachment>gainSlideAttach;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (FdplateReverbPluginAudioProcessorEditor)
};
