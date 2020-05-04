/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
FdplateReverbPluginAudioProcessor::FdplateReverbPluginAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", AudioChannelSet::stereo(), true)
                     #endif
                       ),
#endif
parameters(*this, nullptr, "ParamTreeExample",
{
    std::make_unique<AudioParameterFloat>("gain", "Gain", NormalisableRange<float> (0.0f, 1.0f), 0.5f)
})
{
    gainParam = parameters.getRawParameterValue("gain");
}

FdplateReverbPluginAudioProcessor::~FdplateReverbPluginAudioProcessor()
{
}

//==============================================================================
const String FdplateReverbPluginAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool FdplateReverbPluginAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool FdplateReverbPluginAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool FdplateReverbPluginAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double FdplateReverbPluginAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int FdplateReverbPluginAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int FdplateReverbPluginAudioProcessor::getCurrentProgram()
{
    return 0;
}

void FdplateReverbPluginAudioProcessor::setCurrentProgram (int index)
{
}

const String FdplateReverbPluginAudioProcessor::getProgramName (int index)
{
    return {};
}

void FdplateReverbPluginAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void FdplateReverbPluginAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
}

void FdplateReverbPluginAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool FdplateReverbPluginAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void FdplateReverbPluginAudioProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    
    // In case we have more outputs than inputs, this code clears any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    // This is here to avoid people getting screaming feedback
    // when they first compile a plugin, but obviously you don't need to keep
    // this code if your algorithm always overwrites all the output channels.
    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());

    // This is the place where you'd normally do the guts of your plugin's
    // audio processing...
    // Make sure to reset the state if your inner loop is processing
    // the samples and the outer loop is handling the channels.
    // Alternatively, you can process the samples with the channels
    // interleaved by keeping the same state.
//    buffer.applyGain(*gain);
//    buffer.applyGain(*gainParam);
    buffer.applyGain(*parameters.getRawParameterValue("gain"));
    std::cout << *parameters.getRawParameterValue("gain") << '\n';
        for (int channel = 0; channel < totalNumInputChannels; ++channel)
    {
        auto* channelData = buffer.getWritePointer (channel);

        // ..do something to the data...
    }
}

//==============================================================================
bool FdplateReverbPluginAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* FdplateReverbPluginAudioProcessor::createEditor()
{
    return new FdplateReverbPluginAudioProcessorEditor (*this);
}

//==============================================================================
void FdplateReverbPluginAudioProcessor::getStateInformation (MemoryBlock& destData)
{
//     MemoryOutputStream (destData, true).writeFloat (*gain);
}

void FdplateReverbPluginAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
//    *gain = MemoryInputStream (data, static_cast<size_t> (sizeInBytes), false).readFloat();
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new FdplateReverbPluginAudioProcessor();
}
